#!/usr/bin/env python

import os
import sys
import logging
import pychemia
from pychemia.external.pymatgen import *
from pychemia import pcm_log


def help_info():
    print(""" Create or Update a Database from CIFs

   Use:

       CreateDatabaseFromCIFs.py --dbname 'MongoDB Database name' --path 'Directory with CIFs'
                                [--host localhost ] [--port 27017] [--ssl]
                                [--user None ] [--passwd None]

   """)


version = 0.1

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
            elif option == 'path':
                path = sys.argv[i + 1]
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

    cifs = [x for x in os.listdir(path) if x[-3:] == 'cif']

    for i in cifs:
        pcm_log.info('Reading CIF         : %s' % i)
        structs = cif2structure(path + os.sep + i)
        pcm_log.info('Number of structures: %d' % len(structs))
        for j in structs:
            if j.is_perfect:
                pcm_log.info('Composition         : %s' % str(j.composition))
                pcdb.insert(structure=j)
            else:
                pcm_log.info('DISCARDED: Structure is not perfect')
