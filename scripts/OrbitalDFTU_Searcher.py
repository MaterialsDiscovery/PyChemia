#!/usr/bin/env python

import os
import argparse
import logging
import pychemia


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Orbital DFT+U Searcher Manager')

    parser.add_argument('-host', type=str, help='Hostname or IP of the Mongo Server (default: localhost)',
                        required=False, metavar='hostname')
    parser.add_argument('-dbname', type=str, help='Name of Database', required=True, metavar='name')
    parser.add_argument('-ssl', help='Use SSL (default: False)', action='store_true')
    parser.add_argument('-generation_size', type=int, help='Generation Size (default: 32)', metavar='N', default=32)
    parser.add_argument('-abinit_input', type=str, help='Path to Abinit input file', metavar='path',
                        default='abinit.in')
    parser.add_argument('-new', help='Create new database (default: False)', action='store_true')
    parser.add_argument('-debug', help='Activate debug mode (default: False)', action='store_true')
    parser.add_argument('-clean', help='Clean database before start (default: False)', action='store_true')

    args = parser.parse_args()

    if args.debug:
        loglevel = logging.DEBUG
    else:
        loglevel = logging.INFO

    logging.basicConfig(level=loglevel)
    logger = logging.getLogger('pychemia')
    logger.addHandler(logging.NullHandler())
    logger.setLevel(loglevel)

    if args.host is None:
        host = 'localhost'
    else:
        host = args.host

    print("Host:                  %s" % host)
    print("Database Name:         %s" % args.dbname)
    print("SSL:                   %s" % args.ssl)
    print("Create new database:   %s" % args.new)
    print("Generation Size:       %s" % args.generation_size)
    print("Abinit input:          %s" % args.abinit_input)

    if not os.path.exists(args.abinit_input):
        raise ValueError("ERROR: File %s not found" % args.abinit_input)

    if args.new:
        print("Creating new database")
        admin_name = input('Admin Username:')
        admin_passwd = input('Admin Password:')

    user = input('Username:')
    if user == '':
        user = None
        passwd = None
    else:
        passwd = input('Password:')

    if args.ssl is None:
        ssl = False
    else:
        ssl = args.ssl

    if args.new:
        pcdb = pychemia.db.create_database(name=args.dbname, admin_name=admin_name, admin_passwd=admin_passwd,
                                           user_name=user, user_passwd=passwd, host=host, ssl=args.ssl)
    else:
        dbsettings = {'host': host, 'name': args.dbname, 'user': user, 'passwd': passwd}
        pcdb = pychemia.db.get_database(dbsettings)

    if args.clean:
        pcdb.clean()

    popu = pychemia.population.orbitaldftu.OrbitalDFTU(pcdb, input_path=args.abinit_input)
    searcher = pychemia.searcher.FireFly(popu, generation_size=args.generation_size)
    print("Starting search...")
    searcher.run()
