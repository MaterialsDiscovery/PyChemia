#!/usr/bin/env python
import sys
import logging
import itertools
import argparse
from multiprocessing import Pool, cpu_count
import pychemia
from pychemia.utils.periodic import atomic_symbol

try:
    from qmpy import Entry
except ImportError:
    Entry = None
    print("Could not import 'qmpy' as needed to interface with the OQMD database")
    exit(1)


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

    if 'standard' in a.calculations:
        best_calculation = a.calculations['standard']
    else:
        print('No standard calculation for %3d Calculations available are: %s' % (a.id, a.calculations.keys()))

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

    symm = pychemia.crystal.CrystalSymmetry(structure)
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
    print('Process: %2d Processing from %6d to %6d total: %d' % (start, initial, final, len(entry_ids)))

    for a_id in entry_ids[initial:final]:

        if a_id not in current[index:]:
            ret.append(a_id)
            n += 1
        else:
            index = current.index(a_id)
            # Removing duplicated entries
            if index + 1 < len(current)  and current[index + 1] == a_id:
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
        structure = None
        properties = None
        a = Entry.objects.get(id=entry_id)
        if a.calculation_set.count() > 0:
            structure, properties = run_one(a)
        if structure is not None:
            pcdb.insert(structure, properties=properties)


def getter_star(a_b):
    return getter(*a_b)

version = 0.1
jump = 10000


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Create or Update a PyChemia Database from the OQMD Database (www.oqmd.org)')
    parser.add_argument('-dbname', metavar='<DATABASE>', type=str, help='Database Name', default='PyChemiaDB')
    parser.add_argument('-port', metavar='<PORTNUMBER>', type=int, help='Port (default: 27017)', default=27017)
    parser.add_argument('-ssl', metavar='<SSL>', type=bool, help='Using SSL (default:no)', default=False)
    parser.add_argument('-user', metavar='<USERNAME>', type=str, help='Database Username', default=None)
    parser.add_argument('-host', metavar='<HOSTNAME>', type=str, help='Hostname (default: localhost)', default='localhost')
    parser.add_argument('-nprocs', metavar='N', type=str, help='Number of concurrent proccess (default: Number of CPUs)', default=None)

    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger('pychemia')
    logger.addHandler(logging.NullHandler())
    logger.setLevel(logging.INFO)

    db_settings = {'name': args.dbname, 'host': args.host, 'port': args.port, 'ssl': args.ssl}
    if args.user is not None:
        passwd = input('Password:')
        db_settings['user'] = args.user
        db_settings['passwd'] = passwd

    print('Database settings: \n')
    print(db_settings)
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

    if args.nprocs is None:
        nprocs = cpu_count()
    else:
        nprocs = args.nprocs
    print('Creating a pool of %d processes for feeding the database' % nprocs)
    pool = Pool(processes=nprocs)

    argus = []
    a_args = range((len(entry_ids) / jump) + 1)

    # ret = pool.map(getter_star, itertools.izip(itertools.repeat(entry_ids),
    #                                            itertools.repeat(db_settings),
    #                                            itertools.repeat(current), a_args), chunksize=1)
    #
    # to_insert = []
    # for i in ret:
    #     for j in i:
    #         to_insert.append(j)
    #
    # pool.close()
    #
    # setter(pcdb, to_insert)
