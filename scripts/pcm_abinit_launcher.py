#!/usr/bin/env python

import logging
import os
import sys
import json
import time
import getopt
import pychemia
from pychemia.utils.computing import get_int

logging.basicConfig(level=logging.DEBUG)

try:
    from pychemia.symm import symmetrize
except ImportError:
    def symmetrize(x):
        return x


def cleaner():
    files = ['abinit_LOG', 'abinit_STATUS']
    for j in files:
        for i in os.listdir('.'):
            if i.startswith(j):
                os.remove(i)
    if os.path.exists('abinit.err') and os.path.getsize('abinit.err') == 0:
        os.remove('abinit.err')


def usage(name):
    print("""
NAME
    %s

DESCRIPTION
    ABINIT Run with automatic selection of Norm-Conserving PSPs

OPTIONS

    --help, -h
        Return information on the options and use of this script

    --input, -i <string> (Default: 'abinit.in')
        Structure for which the strength will be computed

    --output, -o <string> (Default: 'pcm_results.json')
        File in JSON format where the results are stored

    --nparal, -n <int> (Default: 2)
        Number of MPI parallel processes for ABINIT

    --binary, -b <string> (Default: 'abinit')
        Path to the ABINIT executable

    --pseudo, -p <string> (Default: 'FHI_LDA')
        Kind and Exchange for Pseudopotentials

""" % os.path.basename(name))


def main(argv):
    try:
        opts, args = getopt.getopt(argv[1:], "hi:o:n:b:p:", ["help", "input=", "output=", "nparal=", "binary=",
                                                             'pseudo='])
    except getopt.GetoptError:
        usage(argv[0])
        sys.exit(2)

    if len(opts) == 0:
        usage(argv[0])
        sys.exit(2)

    # Default Values
    workdir = '.'
    inpu_file = 'abinit.in'
    output_file = 'pcm_results.json'
    nparal = 2
    binary = 'abinit'
    pseudo = 'FHI_LDA'

    for opt, arg in opts:

        if opt in ("-h", "--help"):
            usage(argv[0])
            sys.exit()
        elif opt in ("-i", "--input"):
            inpu_file = arg
        elif opt in ("-o", "--output"):
            output_file = arg
        elif opt in ("-n", "--nparal"):
            nparal = get_int(arg)
        elif opt in ("-b", "--binary"):
            binary = arg
        elif opt in ("-p", "--pseudo"):
            pseudo = arg

    structure = pychemia.code.abinit.AbinitInput(inpu_file).get_structure()
    if structure is None:
        print(" ERROR: Invalid structure, no structure could be obtained")
        sys.exit(2)

    print("\n PyChemia ABINIT Job")
    print(" =======================\n")
    print(" MPI number of procs : ", nparal)
    print(" ABINIT binary       : ", binary)
    print(" ABINIT input_file   : ", inpu_file)
    print(" ABINIT pseudo       : ", pseudo)
    print(" Structure           :\n")
    print(structure)

    wf = open(output_file, 'w')
    data = {'input': {'binary': binary, 'nparal': nparal, 'structure': structure.to_dict}, 'complete': False, }
    json.dump(data, wf, sort_keys=True)
    wf.close()

    cleaner()

    aj = pychemia.code.abinit.AbinitJob()
    aj.initialize(workdir=workdir, input_file='abinit.in',
                  psp_kind=pseudo.split('_')[0],
                  psp_exchange=pseudo.split('_')[1])
    aj.set_inputs()

    if os.path.exists('abinit.out'):
        os.remove('abinit.out')

    aj.run(use_mpi=True, omp_max_threads=2, mpi_num_procs=nparal)
    time.sleep(5)
    ao = pychemia.code.abinit.AbinitOutput('abinit.out')
    old_niter = 0
    while aj.runner.returncode is None:
        aj.runner.poll()
        ao.reload()
        energetics = ao.get_energetics()
        if energetics is not None:
            if len(energetics['iter']) > old_niter:
                print("ITER %4d  ETOT= %15.8f deltaE(h)= %9.2E residm= %9.2E nres2= %9.2E" % (energetics['iter'][-1],
                                                                                              energetics['etot'][-1],
                                                                                              energetics['deltaEh'][-1],
                                                                                              energetics['residm'][-1],
                                                                                              energetics['nres2'][-1]))
                old_niter = len(energetics['iter'])
        time.sleep(60)

    ao.reload()
    energetics = ao.get_energetics()
    occupation_matrix = ao.get_occupation_matrix()
    output = {'energetics': energetics}
    if occupation_matrix is not None:
        output['occupation_matrix'] = occupation_matrix
    if os.path.isfile('abinit-o_OUT.nc'):
        nco = pychemia.code.abinit.netcdf2dict('abinit-o_OUT.nc')
        output['variables'] = nco
    data['output'] = output

    data['complete'] = True

    wf = open(output_file, 'w')
    json.dump(data, wf, sort_keys=True)
    wf.close()

    cleaner()


if __name__ == "__main__":
    main(sys.argv)
