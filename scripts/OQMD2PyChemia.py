#!/usr/bin/env python

import os
import sys
import time
import logging
import itertools
import argparse
from multiprocessing import Pool, cpu_count, Process
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

    energy = 1E10
    best_calculation = None

    if a.calculation_set.count() > 0:

        if 'standard' in a.calculations:
            best_calculation = a.calculations['standard']
            calculation_name = 'standard'
        elif 'fine_relax' in a.calculations:
            best_calculation = a.calculations['fine_relax']
            calculation_name = 'fine_relax'
        elif 'coarse_relax' in a.calculations:
            best_calculation = a.calculations['coarse_relax']
            calculation_name = 'coarse_relax'
        elif 'static' in a.calculations:
            best_calculation = a.calculations['static']
            calculation_name = 'static'
        elif 'relaxation' in a.calculations:
            best_calculation = a.calculations['relaxation']
            calculation_name = 'relaxation'
        elif len(a.calculations) > 0:
            calculations = sorted(a.calculations.keys())
            print('Calculations found: %s, using the last one' % calculations)
            best_calculation = a.calculations[calculations[-1]]
            calculation_name = calculations[-1]
        else:
            print('ERROR: Count > 0 and no calculation found')

    if best_calculation is not None:
        structure_name = None
        if best_calculation.output is not None:
            structure_used = best_calculation.output
            structure_id = best_calculation.output_id
            from_output = True
        elif best_calculation.input is not None:
            print(
                'WARNING: No data was found from the output of the calculation, using input geometries and leaving '
                'energetics empty')
            structure_used = best_calculation.input
            structure_id = best_calculation.input_id
            from_output = False
    else:
        calculation_name = None
        if a.structures is not None and len(a.structures) > 0:
            struct_keys = sorted(a.structures.keys())
            print("WARNING: Calculation not found for %s. Structures found: %s using the first one " % (a, struct_keys))
            structure_used = a.structures[struct_keys[0]]
            structure_id = None
            from_output = False
            structure_name = struct_keys[0]
        else:
            print("ERROR: No calculation and no structure found for %s" % a)
            return None, None

    cell = structure_used.cell.T
    symbols = atomic_symbol(structure_used.atomic_numbers)
    reduced = structure_used.coords

    structure = pychemia.Structure(cell=cell, symbols=symbols, reduced=reduced)
    entry_id = a.id

    if best_calculation is not None:
        calculation_id = best_calculation.id
        energy_pa = best_calculation.energy_pa
        energy = best_calculation.energy
        band_gap = best_calculation.band_gap
        settings = best_calculation.settings

        try:
            spacegroup_number = best_calculation.output.spacegroup.number
        except AttributeError:
            spacegroup_number = None

    else:
        calculation_id = None
        energy_pa = None
        energy = None
        band_gap = None
        settings = None
        spacegroup_number = None
        from_output = False

    try:
        symm = pychemia.crystal.CrystalSymmetry(structure)
        sym2 = symm.number(1E-2)
    except ValueError:
        sym2 = None

    properties = {'oqmd': {'structure_id': structure_id,
                           'entry_id': entry_id,
                           'calculation_id': calculation_id,
                           'energy_pa': energy_pa,
                           'energy': energy,
                           'band_gap': band_gap,
                           'settings': settings,
                           'from_output': from_output,
                           'calculation_name': calculation_name,
                           'structure_name': structure_name,
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
            if index + 1 < len(current) and current[index + 1] == a_id:
                print('We found at least one duplicate!')
                duplicate = False
                for entry in pcdb.db.pychemia_entries.find({'properties.oqmd.entry_id': a_id}):
                    if duplicate:
                        print('Removing PyChemiaDB entry: %s' % str(entry['_id']))
                        pcdb.db.pychemia_entries.remove({'_id': entry['_id']})
                    duplicate = True

    print('Process: %2d Entries missing: %3d' % (start, len(ret)))
    return ret


def setter(db_settings, to_insert):
    print('Processing %d entries - ' % len(to_insert), end='')
    pcdb = pychemia.db.get_database(db_settings)

    if hasattr(os, 'getppid'):  # only available on Unix
        print('parent process: %d - ' % os.getppid(), end='')
    print('process id: %d' % os.getpid())

    index = 0
    for oqmd_id in to_insert:
        if index % 2000 == 0:
            print(index, oqmd_id)
        index += 1
        structure = None
        properties = None

        a = Entry.objects.get(id=oqmd_id)
        structure, properties = run_one(a)

        if structure is not None:
            entry_id = '%d_%s_' % (structure.nspecies, structure.formula)
            n = len(entry_id)
            texto = '%0' + ('%d' % (28 - n)) + 'd'
            entry_id += texto % properties['oqmd']['entry_id']
            if n > 17:
                print("%2d - %s" % (28 - n, entry_id))
            pcdb.insert(structure, properties=properties, entry_id=entry_id)
    return 0


def getter_star(a_b):
    return getter(*a_b)


def setter_star(a_b):
    return setter(*a_b)


version = 0.1
jump = 10000

if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description='Create or Update a PyChemia Database from the OQMD Database (www.oqmd.org)')
    parser.add_argument('-dbname', metavar='<DATABASE>', type=str, help='Database Name', default='PyChemiaMasterDB')
    parser.add_argument('-port', metavar='<PORTNUMBER>', type=int, help='Port (default: 27017)', default=27017)
    parser.add_argument('-ssl', metavar='<SSL>', type=bool, help='Using SSL (default:no)', default=False)
    parser.add_argument('-user', metavar='<USERNAME>', type=str, help='Database Username', default=None)
    parser.add_argument('-host', metavar='<HOSTNAME>', type=str, help='Hostname (default: localhost)',
                        default='localhost')
    parser.add_argument('-nprocs', metavar='N', type=int,
                        help='Number of concurrent proccess (default: Number of CPUs)', default=None)

    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger('pychemia')
    logger.addHandler(logging.NullHandler())
    logger.setLevel(logging.INFO)

    db_settings = {'name': args.dbname, 'host': args.host, 'port': args.port, 'ssl': args.ssl}
    if args.user is not None:
        passwd = input('Password: ')
        db_settings['user'] = args.user
        db_settings['passwd'] = passwd

    print('Database settings: \n%s\b' % db_settings)
    pcdb = pychemia.db.get_database(db_settings)

    nitems = pcdb.entries.count()
    print('Number of entries in the current PyChemia Database: %d' % nitems)

    current = []
    for entry in pcdb.db.pychemia_entries.find({'properties.oqmd.entry_id': {'$exists': True}}):
        current.append(entry['properties']['oqmd']['entry_id'])
    current.sort()
    print('Number of entries coming from OQMD: %d' % len(current))

    print('Number of entries in OQMD...', end='')
    queryset = Entry.objects.all()
    entry_ids = [entry.id for entry in queryset]
    print('%d' % len(entry_ids))

    if args.nprocs is None:
        nprocs = cpu_count()
    else:
        nprocs = args.nprocs

    print('Creating a pool of %d processes for feeding the database' % nprocs)
    pool = Pool(processes=nprocs)

    argus = []
    a_args = range((len(entry_ids) / jump) + 1)

    to_insert = pool.map(getter_star, itertools.izip(itertools.repeat(entry_ids),
                                                     itertools.repeat(db_settings),
                                                     itertools.repeat(current), a_args), chunksize=1)

    pool.close()

    #    to_insert=to_insert[:20]
    print(len(to_insert))
    print(db_settings)

    ps = [None for x in range(nprocs)]
    counter = 0

    # QMPY does not support concurrent executions
    # while counter < len(to_insert):
    #     for i in range(nprocs):
    #         if ps[i] is None or not ps[i].is_alive():
    #             ps[i] = Process(target=setter, args=(db_settings,to_insert[counter]))
    #             ps[i].start()
    #             #ps[i].join()
    #             counter+=1
    #             if counter == len(to_insert):
    #                 break
    #             time.sleep(1)
    #     time.sleep(30)
    # print(ps)

    #    pool = Pool(processes=nprocs)
    #    ret = pool.map(setter_star, itertools.izip(itertools.repeat(db_settings), to_insert), chunksize=1)
    #    pool.close()

    for i in range(len(to_insert)):
        if len(to_insert[i]) > 0:
            setter(db_settings, to_insert[i])
