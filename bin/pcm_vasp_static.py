#!/usr/bin/env python

import logging
logging.basicConfig(level=logging.DEBUG)
import os
import sys
import json
import getopt
import pychemia
from pychemia.code.vasp import ConvergenceKPointGrid, VaspRelaxator, read_poscar, ConvergenceCutOffEnergy
from pychemia.utils.computing import get_float, get_int

try:
    from pychemia.symm import symmetrize
except ImportError:
    def symmetrize(x):
        return x


def cleaner():
    files = ['DOSCAR', 'EIGENVAL', 'IBZKPT', 'OSZICAR', 'PCDAT', 'PROCAR', 'WAVECAR', 'XDATCAR']
    [os.remove(ifile) for ifile in files if os.path.exists(ifile)]


def usage(name):
    print """
NAME
    %s

DESCRIPTION
    Static calculation in VASP with convergence in K-Point Grid and Cut-Off Energy (ENCUT)

OPTIONS

    --help, -h
        Return information on the options and use of this script

    --structure, -s <string> (Default: 'POSCAR')
        Structure for which the strength will be computed

    --output, -o <string> (Default: 'IdealStrength.json')
        File in JSON format where the results are stored

    --nparal, -n <int> (Default: 2)
        Number of MPI parallel processes for VASP

    --binary, -b <string> (Default: 'vasp')
        Path to the VASP executable

    --energy_tol, -e <float> (Default: 1E-3)
        Energy tolerance for three succesive values of the Kpoint grid or CutOff energy

""" % os.path.basename(name)


def main(argv):
    try:
        opts, args = getopt.getopt(argv[1:], "hs:o:n:b:e:", ["help", "structure=", "output=", "nparal=",
                                                             "binary=", "energy_tol="])
    except getopt.GetoptError:
        usage(argv[0])
        sys.exit(2)

    if len(opts) == 0:
        usage(argv[0])
        sys.exit(2)

    # Default Values
    workdir = '.'
    structure_file = 'POSCAR'
    output_file = 'static_calculation.json'
    nparal = 2
    energy_tol = 1E-3
    binary = 'vasp'

    for opt, arg in opts:

        if opt in ("-h", "--help"):
            usage(argv[0])
            sys.exit()
        elif opt in ("-s", "--structure"):
            structure_file = arg
        elif opt in ("-o", "--output"):
            output_file = arg
        elif opt in ("-n", "--nparal"):
            nparal = get_int(arg)
        elif opt in ("-b", "--binary"):
            binary = arg
        elif opt in ("-e", "--energy_tol"):
            energy_tol = get_float(arg)
        elif opt in ("-s", "--structure"):
            structure_file = get_int(arg)

    structure = pychemia.structure_from_file(structure_file)
    if structure is None:
        print " ERROR: Invalid structure, no structure could be obtained from '%s'" % structure_file
        sys.exit(2)

    if structure_file == 'POSCAR':
        os.rename('POSCAR', 'POSCAR_original')

    print "\n PyChemia VASP Static"
    print " =======================\n"
    print " Executable          : ", binary
    print " Energy tolerance    : ", energy_tol
    print " MPI number of procs : ", nparal
    print " Structure           :\n"
    print structure

    wf = open(output_file, 'w')
    data = {'input': {'binary': binary,
                      'energy_tol': energy_tol,
                      'nparal': nparal,
                      'structure': structure.to_dict}}
    json.dumps(data, wf)
    wf.close()

    cleaner()
    print '\nConvergence of Cut-off Energy'
    print '-----------------------------\n'
    ce = ConvergenceCutOffEnergy(structure, energy_tolerance=energy_tol, nparal=nparal)
    ce.run()
    ce.save()
    ce.plot()
    encut = ce.best_encut

    data['output'] = {'ENCUT': encut}
    wf = open(output_file, 'w')
    json.dumps(data, wf)
    wf.close()

    cleaner()
    print '\nConvergence of K-Point Grid'
    print '---------------------------\n'
    ck = ConvergenceKPointGrid(structure, encut=encut, nparal=nparal, energy_tolerance=energy_tol)
    ck.run()
    ck.save()
    ck.plot()
    kp = ck.best_kpoints

    data['output'] = {'KPOINTS': list(kp.grid)}
    wf = open(output_file, 'w')
    json.dumps(data, wf)
    wf.close()

    cleaner()

if __name__ == "__main__":

    main(sys.argv)
