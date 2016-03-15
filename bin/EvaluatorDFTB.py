#!/usr/bin/env python

import logging
import argparse
from pychemia.code.dftb import EvaluatorDaemon


if __name__ == '__main__':

    logging.basicConfig(level=logging.DEBUG)
    logger = logging.getLogger('pychemia')
    logger.addHandler(logging.NullHandler())
    logger.setLevel(logging.DEBUG)

    description = """Launch DFTB+ for non-evaluated entries in a PyChemia Database"""

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
                        default=None, metavar='dbname', type=str,
                        help='PyChemia Database name (default: None)')
    parser.add_argument('-l', '--slater_path',
                        default=None, metavar='path', type=str,
                        help='Slater Path (default: None)')
    parser.add_argument('-f', '--target_forces',
                        default=1E-3, metavar='x', type=float,
                        help='Target Forces (default: 1E-3)')
    parser.add_argument('-n', '--nparal',
                        default=1, metavar='N', type=int,
                        help='Number of parallel processes (default: 1)')
    parser.add_argument('-r', '--replicaset',
                        default=None, metavar='name', type=str,
                        help='ReplicaSet  (default: None)')
    parser.add_argument('-w', '--workdir',
                        default='.', metavar='path', type=str,
                        help='Working Directory  (default: None)')
    parser.add_argument('--evaluate_all', action='store_true',
                        help='Evaluate All  (default: No)')
    parser.add_argument('--waiting', action='store_true',
                        help='Waiting  (default: No)')
    parser.add_argument('--ssl', action='store_true',
                        help='Use SSL to connect to MongoDB  (default: No)')

    args = parser.parse_args()

    if args.dbname is None:
        parser.print_help()
        exit(1)

    print args

    db_settings = {'name': args.dbname, 'host': args.host, 'port': args.port, 'ssl': args.ssl,
                   'replicaset': args.replicaset}
    if args.user is not None:
        if args.passwd is None:
            raise ValueError('Password is mandatory if user is entered')
        db_settings['user'] = args.user
        db_settings['passwd'] = args.passwd
    relaxator_params = {'slater_path': args.slater_path}
    print 'pyChemia Evaluator using DFTB+'
    print 'dbname    : %s' % args.dbname
    print 'host      : %s' % args.host
    print 'port      : %d' % args.port
    print 'user      : %s' % args.user
    print 'replicaset: %s' % args.replicaset
    print 'workdir   : %s' % args.workdir
    print 'nparal    : %d' % args.nparal
    print 'slater-path   : %s' % str(args.slater_path)
    print 'target-forces : %.2E' % args.target_forces
    print 'evaluate_all  : %s' % str(args.evaluate_all)
    print 'waiting       : %s' % str(args.waiting)
    print 'ssl           : %s' % str(args.ssl)

    print db_settings

    evaluator = EvaluatorDaemon(db_settings, args.workdir,
                                target_forces=args.target_forces,
                                nparal=args.nparal,
                                relaxator_params=relaxator_params,
                                evaluate_all=args.evaluate_all,
                                waiting=args.waiting)
    evaluator.run()
