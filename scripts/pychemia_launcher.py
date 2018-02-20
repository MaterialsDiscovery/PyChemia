#!/usr/bin/env python

import os
import sys
import logging
import pychemia
import subprocess
import xml.etree.ElementTree as ElementTree

logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger('pychemia')
logger.addHandler(logging.NullHandler())
logger.setLevel(logging.DEBUG)

version = 0.1


def help_info():
    print(""" pychemia_launcher.py
Search for suitable directories for relaxation

   Use:

       pychemia_launcher.py --user user_name
                            --queue queue_name
                            --mail mail_job_info
                            --path where_to_search_for_structures
                            [--binary name_of_executable]
   """)


def find_structures(basedir, filename='to_relax'):
    ret = []
    if filename in os.listdir(basedir):
        ret.append(basedir)
    for i in [x for x in os.listdir(basedir) if os.path.isdir(basedir + os.sep + x)]:
        ret += find_structures(basedir + os.sep + i)
    return ret


def get_structure_file(basedir):
    files = os.listdir(basedir)
    for ifile in ['structure_phase2.json', 'structure_phase1.json', 'structure2.json',
                  'structure1.json', 'structure.json', 'POSCAR']:
        if ifile in files:
            return ifile
    return None


if __name__ == '__main__':

    # Script starts from here
    if len(sys.argv) < 2:
        help_info()
        sys.exit(1)

    user = None
    mail = None
    queue = None
    path = None
    binary = 'vasp'
    nparal = 4
    nhours = 24

    for istruct in range(1, len(sys.argv)):
        if sys.argv[istruct].startswith('--'):
            option = sys.argv[istruct][2:]
            # fetch sys.argv[1] but without the first two characters
            if option == 'version':
                print(version)
                sys.exit()
            elif option == 'help':
                help_info()
                sys.exit()
            elif option == 'user':
                user = sys.argv[istruct + 1]
            elif option == 'mail':
                mail = sys.argv[istruct + 1]
            elif option == 'queue':
                queue = sys.argv[istruct + 1]
            elif option == 'path':
                path = sys.argv[istruct + 1]
            elif option == 'binary':
                binary = sys.argv[istruct + 1]
            elif option == 'nparal':
                nparal = int(sys.argv[istruct + 1])
            elif option == 'nhours':
                nhours = int(sys.argv[istruct + 1])
            else:
                print('Unknown option. --' + option)

    if path is None or user is None or queue is None:
        help_info()

    data = subprocess.check_output(['qstat', '-x', '-f', '-u', user])
    xmldata = ElementTree.fromstring(data)
    jobs = [i.find('Job_Name').text for i in xmldata.findall('Job')]

    found_structures = find_structures(path)

    for istruct in found_structures:
        if os.path.isfile(istruct + os.sep + 'lock'):
            print("Locked:    %s" % istruct)
        elif os.path.basename(istruct) in jobs:
            print('Submitted: %s' % istruct)
        else:
            print('To submit: %s' % istruct)
            rr = pychemia.runner.PBSRunner(istruct)
            rr.initialize(nodes=1, ppn=nparal, mail=mail, queue=queue, walltime=[0, nhours, 0, 0])
            structure_file = get_structure_file(istruct)
            if structure_file is None:
                print('No suitable structure was found on: ', istruct)
            rr.set_template('pcm_vasp_relaxator.py --binary %s --nparal %d --structure_file %s'
                            % (binary, nparal, structure_file))
            rr.write_pbs()
            rr.submit()
