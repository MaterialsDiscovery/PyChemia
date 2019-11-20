#!/usr/bin/env python

import os
import logging
import argparse
import time
import pychemia
from pychemia.runner import get_jobs

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
    print('')
    print('dbname    : %s' % args.dbname)
    print('source_dir: %s' % args.source_dir)
    print('')
    print('host      : %s' % args.host)
    print('port      : %d' % args.port)
    print('user      : %s' % args.user)
    print('replicaset: %s' % args.replicaset)
    print('binary    : %s' % str(args.binary))
    print('ssl       : %s' % str(args.ssl))

    assert(len(args.dbname) == len(args.source_dir))

    while True:

        for idb in range(len(args.dbname)):

            print('DATABASE: %s' % args.dbname[idb])
            db_settings['name'] = args.dbname[idb]
            pcdb = pychemia.db.get_database(db_settings)
            popu = pychemia.population.NonCollinearMagMoms(pcdb, source_dir=args.source_dir[idb])

            print('Number of candidates evaluated: %d' % len(popu.actives_evaluated))
            print('Number of candidates not evaluated: %d' % len(popu.actives_no_evaluated))

            to_compute = popu.actives_no_evaluated

            print('Candidates to compute:')
            for ijob in to_compute:
                print(ijob)
            jobs = get_jobs(args.pbs_user)
            jobnames = [jobs[x]['Job_Name'] for x in jobs]
            source_dir = args.source_dir[idb]

            for ijob in to_compute:
                if str(ijob) not in jobnames:
                    data_collected = popu.collect_data(ijob, workdir=source_dir + os.sep + str(ijob))
                    if not data_collected:
                        print('Preparing and submitting job: %s' % str(ijob))
                        popu.prepare_folder(ijob, workdir=source_dir + os.sep + str(ijob))

                        pbs = pychemia.runner.PBSRunner(source_dir + os.sep + str(ijob))
                        pbs.initialize(ppn=args.pbs_ppn, walltime=[args.pbs_nhours, 0, 0], mail=args.pbs_mail,
                                       queue=args.pbs_queue)
                        pbs.set_template(source_dir+os.sep+'template.pbs')
                        pbs.write_pbs()
                        pbs.submit()
                else:
                    print('Job %s is on queue or running' % str(ijob))
        print('I will be waiting for 60 minutes before checking new candidates')
        time.sleep(3600)
