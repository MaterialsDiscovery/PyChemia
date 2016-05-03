#!/usr/bin/env python
import sys
import logging
import itertools
from multiprocessing import Pool
import pychemia
from pychemia.utils.periodic import atomic_symbol

try:
    from qmpy import Entry
except ImportError:
    Entry = None
    print("Could not import 'qmpy' as needed to interface with the OQMD database")
    exit(1)

version = 0.1
jump = 10000


def help_info():
    print(""" Create or Update a PyChemia Database from the OQMD Database (www.oqmd.org)

   Use:

       OQMD2PyChemia.py --dbname 'MongoDB Database name'
                                [--host localhost ] [--port 27017]
                                [--user None ] [--passwd None] [--ssl]

   """)


def run_one(a):
    """
    Take one OQMD 'Entry' object, search all the calculations associated and take the best calculation
    in order to insert its data into the PyChemia Database

    :param a: OQMD Entry object
    :return:
    """
    print('Entry: %6d  Number of Calculations: %3d  Energies: %s' % (a.id, a.calculation_set.count(),
                                                                     str([c.energy for c in a.calculation_set.all()])))
    energy = 1E10
    best_calculation = None
    for c in a.calculation_set.all():
        if c.energy is not None and c.energy < energy:
            best_calculation = c
            energy = c.energy

    if best_calculation is None:
        return None, None

    cell = best_calculation.output.cell.T
    symbols = atomic_symbol(best_calculation.output.atomic_numbers)
    reduced = best_calculation.output.coords

    structure = pychemia.Structure(cell=cell, symbols=symbols, reduced=reduced)
    structure_id = best_calculation.output_id
    entry_id = a.id
    calculation_id = best_calculation.id
    energy_pa = best_calculation.energy_pa
    energy = best_calculation.energy
    band_gap = best_calculation.band_gap
    settings = best_calculation.settings
    try:
        spacegroup_number = best_calculation.output.spacegroup.number
    except ValueError:
        spacegroup_number = None

    symm = pychemia.symm.StructureSymmetry(structure)
    sym2 = symm.number(1E-2)

    properties = {'oqmd': {'structure_id': structure_id,
                           'entry_id': entry_id,
                           'calculation_id': calculation_id,
                           'energy_pa': energy_pa,
                           'energy': energy,
                           'band_gap': band_gap,
                           'settings': settings,
                           'spacegroup_number': spacegroup_number},
                  'spacegroup_number': {'value': sym2, 'symprec': 1E-2}}

    return structure, properties


def getter(entry_ids, db_settings, current, start=0):
    pcdb = pychemia.db.get_database(db_settings)

    ret = []

    index = 0
    n = 0

    initial = start * jump
    final = min(start * jump + jump, len(entry_ids))
    print('Process: %2d Processing from %6d to %6d' % (start, initial, final))

    for a_id in entry_ids[initial:final]:

        if a_id not in current[index:]:
            ret.append(a_id)
            n += 1
        else:
            index = current.index(a_id)
            # Removing duplicated entries
            if index < len(current) and current[index + 1] == a_id:
                print('We found at least one duplicate!')
                duplicate = False
                for entry in pcdb.db.pychemia_entries.find({'properties.oqmd.entry_id': a_id}):
                    if duplicate:
                        print('Removing PyChemiaDB entry: %s' % str(entry['_id']))
                        pcdb.db.pychemia_entries.remove({'_id': entry['_id']})
                    duplicate = True

    print('Process: %2d Entries missing: %3d' % (start, len(ret)))

    return ret


def setter(pcdb, to_insert):
    for entry_id in to_insert:
        a = Entry.objects.get(id=entry_id)
        structure, properties = run_one(a)
        if structure is not None:
            pcdb.insert(structure, properties=properties)


def getter_star(a_b):
    return getter(*a_b)


if __name__ == '__main__':

    logging.basicConfig(level=logging.DEBUG)
    logger = logging.getLogger('pychemia')
    logger.addHandler(logging.NullHandler())
    logger.setLevel(logging.DEBUG)

    # Script starts from here
    if len(sys.argv) < 2:
        help_info()
        sys.exit(1)

    dbname = None
    host = 'localhost'
    port = 27017
    user = None
    passwd = None
    path = None
    ssl = False

    for i in range(1, len(sys.argv)):
        if sys.argv[i].startswith('--'):
            option = sys.argv[i][2:]
            # fetch sys.argv[1] but without the first two characters
            if option == 'version':
                print(version)
                sys.exit()
            elif option == 'help':
                help_info()
                sys.exit()
            elif option == 'dbname':
                dbname = sys.argv[i + 1]
            elif option == 'host':
                host = sys.argv[i + 1]
            elif option == 'port':
                dbname = int(sys.argv[i + 1])
            elif option == 'user':
                user = sys.argv[i + 1]
            elif option == 'passwd':
                passwd = sys.argv[i + 1]
            elif option == 'ssl':
                ssl = True
            else:
                print('Unknown option. --' + option)

    if dbname is None:
        help_info()
        sys.exit(1)

    db_settings = {'name': dbname, 'host': host, 'port': port, 'ssl': ssl}
    if user is not None:
        if passwd is None:
            raise ValueError('Password is mandatory if user is entered')
        db_settings['user'] = user
        db_settings['passwd'] = passwd

    pcdb = pychemia.db.get_database(db_settings)

    nitems = pcdb.entries.count()
    print('Number of entries in the current PyChemia Database: %d' % nitems)

    current = []
    for entry in pcdb.db.pychemia_entries.find({'properties.oqmd.entry_id': {'$exists': True}}):
        current.append(entry['properties']['oqmd']['entry_id'])
    current.sort()
    print('Number of entries coming from OQMD: %d' % len(current))

    queryset = Entry.objects.all()
    entry_ids = [entry.id for entry in queryset]
    print('Number of entries in OQMD: %d' % len(entry_ids))

    pool = Pool(processes=6)

    argus = []
    a_args = range(len(entry_ids) / jump)

    ret = pool.map(getter_star, itertools.izip(itertools.repeat(entry_ids),
                                               itertools.repeat(db_settings),
                                               itertools.repeat(current), a_args))

    to_insert = []
    for i in ret:
        for j in i:
            to_insert.append(j)

    pool.close()

    setter(pcdb, to_insert)
