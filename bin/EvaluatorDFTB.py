#!/usr/bin/env python

import os
import sys
import logging
from pychemia.code.dftb import EvaluatorDaemon

version = 0.1


def help_info():
    print(""" Evaluator for DFTB+
Start the execution of the Evaluator using DFTB+ as relaxator

   Use:

       EvaluatorDFTB.py --dbname 'MongoDB Database name'
                        [--host localhost ] [--port 27017]
                        [--user None ] [--passwd None] [--ssl]
                        [--slater-path $HOME] [--workdir $HOME]
                        [--target-forces 1E-3] [--nparal 2]
                        [--evaluate-all] [--waiting]
   """)


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
    slater_path = [os.getenv('HOME')]
    workdir = os.getenv('HOME')
    target_forces = 1E-3
    nparal = 2
    evaluate_all = False
    waiting = False
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
            elif option == 'workdir':
                workdir = sys.argv[i + 1]
            elif option == 'slater-path':
                slater_path.append(sys.argv[i + 1])
            elif option == 'nparal':
                nparal = int(sys.argv[i + 1])
            elif option == 'target-forces':
                target_forces = float(sys.argv[i + 1])
            elif option == 'evaluate-all':
                evaluate_all = True
            elif option == 'waiting':
                waiting = True
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
    relaxator_params = {'slater_path': slater_path}
    print 'pyChemia Evaluator using DFTB+'
    print 'dbname  : %s' % dbname
    print 'host    : %s' % host
    print 'port    : %d' % port
    print 'user    : %s' % user
    print 'workdir : %s' % workdir
    print 'nparal  : %d' % nparal
    print 'slater-path   : %s' % str(slater_path)
    print 'target-forces : %.2E' % target_forces
    print 'evaluate_all  : %s' % str(evaluate_all)
    print 'waiting       : %s' % str(waiting)
    print 'ssl           : %s' % str(ssl)

    evaluator = EvaluatorDaemon(db_settings, workdir,
                                target_forces=target_forces,
                                nparal=nparal,
                                relaxator_params=relaxator_params,
                                evaluate_all=evaluate_all,
                                waiting=waiting)
    evaluator.run()
