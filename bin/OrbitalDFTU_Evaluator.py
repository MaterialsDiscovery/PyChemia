#!/usr/bin/env python

from builtins import input

import os
import sys
import argparse
import pychemia


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Orbital DFTU Evaluator and Analysis Tool')

    parser.add_argument('-host', type=str, help='Hostname or IP of the Mongo Server', required=True)
    parser.add_argument('-name', type=str, help='Name of Database', required=True)
    parser.add_argument('-dbuser', type=str, help='Username on Database', required=True)
    parser.add_argument('-quser', type=str, help='Username on PBS Queue ', required=True)
    parser.add_argument('-queue', type=str, help='Queue', required=True)
    parser.add_argument('-hours', type=int, help='Number of hours to request', required=True)
    parser.add_argument('-features', type=str, help='Features for queue', required=False, default=None)
    parser.add_argument('-basepath', type=str, help='Path where calculations are performed', default='.')

    args = parser.parse_args()
    basepath = args.basepath    

    passwd = input('Password:')

    if not os.path.isdir(args.basepath) or not os.path.isfile(basepath+'/abinit.in'):
        print('ERROR: Wrong basepath %s' % basepath)
        parser.print_help()
        sys.exit(1)

    db_settings = {'name': args.name, 'host': args.host, 'user': args.dbuser, 'passwd': passwd}
    pcdb = pychemia.db.get_database(db_settings)
    popu = pychemia.population.OrbitalDFTU(pcdb, basepath+'/abinit.in')
    popu.evaluator(username=args.quser, basedir=basepath, queue=args.queue,
                   walltime=[args.hours, 0, 0], ppn=6, features=args.features)
