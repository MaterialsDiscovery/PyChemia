#!/usr/bin/env python

import logging
import argparse
import time
import pychemia
from pychemia.runner.pbs import get_jobs

if __name__ == '__main__':

    logging.basicConfig(level=logging.DEBUG)
    logger = logging.getLogger('pychemia')
    logger.addHandler(logging.NullHandler())
    logger.setLevel(logging.DEBUG)

    description = """Launch VASP for non-evaluated entries in a PyChemia Database"""

    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-t', '--host',
                        default='localhost', metavar='server', type=str,
                        help='Hostname or address (default: localhost)')
    parser.add_argument('-o', '--port',
                        default=27017, metavar='port', type=int,
                        help='MongoDB port (default: 27017)')
    parser.add_argument('-u', '--user',
                        default=None, metavar='username', type=str,
                        help='Username (default: None)')
    parser.add_argument('-p', '--passwd',
                        default=None, metavar='password', type=str,
                        help='Password (default: None)')
    parser.add_argument('-d', '--dbname',
                        default=None, metavar='dbname', type=str, nargs='+',
                        help='PyChemia Database name (default: None)')
    parser.add_argument('-b', '--binary',
                        default='vasp', metavar='path', type=str,
                        help='VASP binary (default: vasp)')
    parser.add_argument('-s', '--source_dir',
                        default=None, metavar='path', type=str, nargs='+',
                        help='Source Directory where KPOINTS, POSCAR, INCAR and POTCAR are (default: None)')
    parser.add_argument('-r', '--replicaset',
                        default=None, metavar='name', type=str,
                        help='ReplicaSet  (default: None)')
    parser.add_argument('--ssl', action='store_true',
                        help='Use SSL to connect to MongoDB  (default: No)')
    parser.add_argument('--pbs_ppn',
                        default=1, metavar='N', type=int,
                        help='Number of MPI parallel processes (default: 1)')
    parser.add_argument('--pbs_mail',
                        default=None, metavar='user@mail.server', type=str,
                        help='Mail address for PBS (default: None)')
    parser.add_argument('--pbs_queue',
                        default=None, metavar='queue_name', type=str,
                        help='Queue for PBS (default: None)')
    parser.add_argument('--pbs_nhours',
                        default=24, metavar='N', type=int,
                        help='Number of hours for PBS (default: 24)')
    parser.add_argument('--pbs_user',
                        default=None, metavar='username', type=str,
                        help='Username for PBS (default: None)')

    args = parser.parse_args()
    if args.dbname is None:
        parser.print_help()
        exit(1)

    print(args)
    db_settings = {'host': args.host, 'port': args.port, 'ssl': args.ssl, 'replicaset': args.replicaset}
    if args.user is not None:
        if args.passwd is None:
            raise ValueError('Password is mandatory if user is entered')
        db_settings['user'] = args.user
        db_settings['passwd'] = args.passwd
    print('pyChemia Evaluator using VASP')
    print('dbname    : %s' % args.dbname)
    print('host      : %s' % args.host)
    print('port      : %d' % args.port)
    print('user      : %s' % args.user)
    print('replicaset: %s' % args.replicaset)
    print('binary    : %s' % str(args.binary))
    print('ssl       : %s' % str(args.ssl))


    while True:

        for iname in args.dbname:

            print('DATABASE: %s' % iname)
            db_settings['name'] = iname
            pcdb = pychemia.db.get_database(db_settings)
            popu = pychemia.population.NonCollinearMagMoms(pcdb, source_dir=args.source_dir)

            print('Number of candidates evaluated: %d' % len(popu.actives_evaluated))

            to_compute = popu.actives_no_evaluated

            print('Candidates to compute:')
            for i in to_compute:
                print(i)
            current_jobs = get_jobs(args.pbs_user)

            for i in to_compute:
                if str(i) not in current_jobs:
                    data_collected = popu.collect_data(i, str(i))
                    if not data_collected:
                        print('Preparing and submitting job: %s' % str(i))
                        popu.prepare_folder(i, workdir=str(i))

                        pbs = pychemia.runner.PBSRunner(str(i))
                        pbs.initialize(ppn=args.pbs_ppn, walltime=[args.pbs_nhours, 0, 0], mail=args.pbs_mail,
                                       queue=args.pbs_queue)
                        pbs.set_template('template.pbs')
                        pbs.write_pbs()
                        pbs.submit()
                else:
                    print('Job %s is on queue or running' % str(i))
        time.sleep(3600)
