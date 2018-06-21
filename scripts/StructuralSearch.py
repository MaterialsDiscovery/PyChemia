#!/usr/bin/env python3

import os
import json
import sys
import time
import argparse
import multiprocessing
import pychemia
from pychemia.db import get_database


def safe_read_json(filename):
    """
    Safely read a given filename, extract and returns the associated dictionary
    """
    if not os.path.exists(filename):
        raise ValueError("ERROR: Could not read file: %s" % filename)
    rf = open(filename)
    try:
        data = json.load(rf)
    except ValueError:
        raise ValueError("ERROR: File is not in proper JSON format: %s" % filename)
    rf.close()
    return data


def searcher(dbsettings, popu_settings, searcher_settings, composition, max_formulas):

    pcdb = get_database(dbsettings)

    if 'target_forces' in popu_settings:
        target_forces = popu_settings['target_forces']
    else:
        target_forces = 1E-3
    if 'value_tol' in popu_settings:
        value_tol = popu_settings['value_tol']
    else:
        value_tol = 1E-2
    if 'distance_tolerance' in popu_settings:
        distance_tolerance = popu_settings['distance_tolerance']
    else:
        distance_tolerance = 0.3

    # Population
    popu = pychemia.population.RelaxStructures(pcdb, composition=composition, tag=composition,
                                               target_forces=target_forces, value_tol=value_tol,
                                               distance_tolerance=distance_tolerance, min_comp_mult=1,
                                               max_comp_mult=max_formulas)

    # Global Searcher
    if 'target_value' in searcher_settings:
        target_value = searcher_settings['target_value']
    else:
        target_value = None

    method = searcher_settings['method']
    if method == 'FireFly':
        globsearcher = pychemia.searcher.FireFly(popu, params=searcher_settings['params'],
                                                 generation_size=searcher_settings['generation_size'],
                                                 stabilization_limit=searcher_settings['stabilization_limit'],
                                                 target_value=target_value)

    globsearcher.run()


def manager(population, code):
    pass


def worker(population, code, candidate):
    pass


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='PyChemia Structural Searcher')

    subparsers = parser.add_subparsers(help='commands', dest='subparser_name')

    # The Searcher uses global searchers, one for each composition 
    searcher_parser = subparsers.add_parser('searcher', help='Create one searcher for each composition')

    searcher_parser.add_argument('-clean', action='store_true',
                                 help='If true, clean the entire database before populate',
                                 required=False, default=False)
    searcher_parser.add_argument('--db_settings', '-d', type=str,
                                 help='Path to JSON file with database settings',
                                 required=True)
    searcher_parser.add_argument('--population_settings', '-p', type=str,
                                 help='Path to JSON file with population settings',
                                 required=True)
    searcher_parser.add_argument('--searcher_settings', '-s', type=str,
                                 help='Path to JSON file with searcher settings',
                                 required=True)
    searcher_parser.add_argument('-compositions', '-c', type=str,
                                 help='Path to JSON file with settings for all compositions to search',
                                 required=True)

    # The Manager search the database for unevaluated candidates, create submission jobs for them and 
    # collect results after execution
    manager_parser = subparsers.add_parser('manager', help='Search a database for candidates to be evaluated')

    manager_parser.add_argument('-db_settings', '-d', type=str,
                                help='Path to JSON file with database settings',
                                required=True)
    manager_parser.add_argument('-population_settings', '-p', type=str,
                                help='Path to JSON file with settings relative to relaxation',
                                required=True)
    manager_parser.add_argument('-basepath', type=str,
                                help='Path where calculations are performed',
                                required=False, default='.')

    # The Worker takes one candidate and execute the selected code to relax the structure
    worker_parser = subparsers.add_parser('worker', help='Takes one candidate and compute a local relaxation on it')

    worker_parser.add_argument('-case', type=str,
                               help='ID of candidate to relax',
                               required=True)
    worker_parser.add_argument('-basepath', type=str,
                               help='Path where calculations are performed',
                               required=False, default='.')

    args = parser.parse_args()
    print(args)

    if args.subparser_name in ['searcher', 'manager']:

        db = safe_read_json(args.db_settings)
        pcdb = get_database(db) 

        print("")
        print("Database Settings")
        print("-----------------\n")
        print(pcdb)
        print("")

        popu_settings = safe_read_json(args.population_settings)
        print("Population Settings")
        print("-------------------\n")
        for i in popu_settings:
            print(" %30s : %s" % (i.ljust(30), popu_settings[i]))
        print("")

    if args.subparser_name == 'searcher':

        searcher_settings = safe_read_json(args.searcher_settings)
        print("Searcher Settings")
        print("-----------------\n")
        for i in searcher_settings:
            print(" %30s : %s" % (i.ljust(30), searcher_settings[i]))
        print("")

        comp_settings = safe_read_json(args.compositions)
        print("Compositions and max number of formulas allowed")
        print("---------------------------------------------\n")
        for i in comp_settings:
            print(" %30s : %s" % (i.ljust(30), comp_settings[i]))
        print("")

        # Check consistency of searcher settings
        if 'method' not in searcher_settings:
            print('Not "method" was declared in searcher settings')
            sys.exit(1)
        if 'generation_size' not in searcher_settings:
            print('Not "generation_size" was declared in searcher settings')
            sys.exit(1)
        if 'stabilization_limit' not in searcher_settings:
            print('Not "stabilization_limit" was declared in searcher settings')
            sys.exit(1)

        sprocs = {}
        for composition in comp_settings:
            max_formulas = comp_settings[composition]
            sprocs[composition] = multiprocessing.Process(target=searcher,
                                                          args=(db,
                                                                popu_settings,
                                                                searcher_settings,
                                                                composition,
                                                                max_formulas))
            sprocs[composition].start()
        
    elif args.subparser_name == 'manager':

        popu = pychemia.population.RelaxStructures(pcdb)
        print(popu)

        while True:
            
            for icase in popu.members:
                if not popu.is_evaluated(icase):
                    ipath = args.basepath + os.sep + str(icase)
                    if not os.path.isdir(ipath):
                        os.mkdir(ipath)

                    st = popu.get_structure(icase) 
                    task = pychemia.code.abinit.task.IonRelaxation(st)

