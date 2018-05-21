#!/usr/bin/env python

import logging
import os
import sys
import json
import getopt
import numpy as np
import pychemia
import pandas
from pychemia.utils.computing import get_int, get_float
from pychemia.code.vasp.task import ConvergenceKPointGrid, ConvergenceCutOffEnergy, IdealStrength

logging.basicConfig(level=logging.DEBUG)


def cleaner():
    files = ['DOSCAR', 'EIGENVAL', 'IBZKPT', 'OSZICAR', 'PCDAT', 'PROCAR', 'WAVECAR', 'XDATCAR']
    [os.remove(ifile) for ifile in files if os.path.exists(ifile)]


def usage(name):
    print("""
NAME

    %s

DESCRIPTION
    Computes the ideal strength of a material, the cell is deformed along the axis especified by 'expansion', the
    structure is enlarge of compress from ini_factor to fin_factor with intervals of delta_factor. A convergence
    of Kpoints and Cutoff energy are applied to preserve an accuracy of 'energy_tol'. The VASP code is used

OPTIONS

    --help, -h
        Return information on the options and use of this script

    --structure, -s <string> (Default: 'POSCAR')
        Structure for which the strength will be computed

    --output, -o <string> (Default: 'IdealStrength.json')
        File in JSON format where the results are stored

    --ini_factor, -i <float> (Default: 0.01)
        Initial factor for which the cell will be enlarge or compress. A value of 0.01 means increase by 1 percent

    --fin_factor, -f <float> (Default: 0.1)
        Final factor for enlarging or compressing the cell

    --delta_factor, -d <float> (Default: 0.01)
        Spacing between values of the factor

    --nparal, -n <int> (Default: 2)
        Number of MPI parallel processes for VASP

    --binary, -b <string> (Default: 'vasp')
        Path to the VASP executable

    --energy_tol, -e <float> (Default: 1E-3)
        Energy tolerance for three succesive values of the Kpoint grid or CutOff energy

    --target_forces, -t <float> (Default: 1E-3)
        Target forces for internal relaxation

    --expansion, -x <int> (Default: 111)
        Determines along which directions the lattice is deformed
        111 means along the 3 directions
        101 means along only along the 'a' and 'c' lattice lengths
""" % os.path.basename(name))


def main(argv):
    try:
        opts, args = getopt.getopt(argv[1:], "hs:o:i:f:d:n:b:e:x:t:", ["help", "structure=", "output=", "ini_factor=",
                                                                       "fin_factor=", "delta_factor", "nparal=",
                                                                       "binary=", "energy_tol=", "expansion=",
                                                                       "target_forces"])
    except getopt.GetoptError:
        usage(argv[0])
        sys.exit(2)

    if len(opts) == 0:
        usage(argv[0])
        sys.exit(2)

    # Default Values
    structure_file = 'POSCAR'
    delta_factor = 0.01
    ini_factor = 0.01
    fin_factor = 0.1
    output_file = 'IdealStrength.json'
    nparal = 2
    energy_tol = 1E-3
    target_forces = 1E-3
    binary = 'vasp'
    expansion = 111

    for opt, arg in opts:
        if opt in ("-h", "--help"):
            usage(argv[0])
            sys.exit()
        elif opt in ("-s", "--structure"):
            structure_file = arg
        elif opt in ("-o", "--output"):
            output_file = arg
        elif opt in ("-i", "--ini_factor"):
            ini_factor = get_float(arg)
        elif opt in ("-f", "--fin_factor"):
            fin_factor = get_float(arg)
        elif opt in ("-d", "--delta_factor"):
            delta_factor = get_float(arg)
        elif opt in ("-n", "--nparal"):
            nparal = get_int(arg)
        elif opt in ("-b", "--binary"):
            binary = arg
        elif opt in ("-e", "--energy_tol"):
            energy_tol = get_float(arg)
        elif opt in ("-t", "--target_forces"):
            target_forces = get_float(arg)
        elif opt in ("-x", "--expansion"):
            expansion = get_int(arg)
        elif opt in ("-s", "--structure"):
            structure_file = get_int(arg)

    expansion = [int(expansion / 100), int(expansion / 10) % 10, expansion % 10]
    if len(expansion) == 0:
        print(" ERROR: ini_factor, fin_factor and delta_factor are not creating a finite range of values")
        sys.exit(2)

    structure = pychemia.structure_from_file(structure_file)
    if structure is None:
        print(" ERROR: Invalid structure, no structure could be obtained from '%s'" % structure_file)
        sys.exit(2)

    if structure_file == 'POSCAR':
        os.rename('POSCAR', 'POSCAR_original')

    print("\n PyChemia Ideal Strenght")
    print(" =======================\n")
    print(" Scaling factors     : ", str(np.arange(ini_factor, fin_factor + 0.9 * delta_factor, delta_factor)))
    print(" Executable          : ", binary)
    print(" Energy tolerance    : ", energy_tol)
    print(" Target forces       : ", target_forces)
    print(" Expansion directions: ", str(expansion))
    print(" MPI number of procs : ", nparal)
    print(" Structure           :\n")
    print(structure)

    cleaner()
    print('\nConvergence of Cut-off Energy')
    print('-----------------------------\n')
    ce = ConvergenceCutOffEnergy(structure, energy_tolerance=energy_tol, binary=binary)
    if os.path.isfile('convergence_encut.json'):
        print('A previous convergence study was found, loading...')
        ce.load()
    if not ce.is_converge:
        ce.run(nparal)
        ce.save()
        ce.plot()
    encut = ce.best_encut
    print('ENCUT: ', encut)

    cleaner()
    print('\nConvergence of K-Point Grid')
    print('---------------------------\n')
    ck = ConvergenceKPointGrid(structure, workdir='.', binary=binary, encut=encut,
                               energy_tolerance=energy_tol, recover=True)
    if os.path.isfile('convergence_kpoints.json'):
        print('A previous convergence study was found, loading...')
        ck.load()
    if not ck.is_converge:
        ck.run(nparal)
        ck.save()
        ck.plot()
    kp = ck.best_kpoints
    kp_density = kp.get_density_of_kpoints(structure.lattice)
    print('KPOINT GRID', kp.grid)

    strenght = IdealStrength(structure, ini_factor, fin_factor, delta_factor, kp, kp_density, expansion, encut,
                             nparal, binary, target_forces, output_file)

    strenght.run(nparal)
    strenght.save()

    cleaner()


def plot_strain(filename):
    import matplotlib.pyplot as plt
    a = json.load(open(filename))
    df = pandas.DataFrame(a)

    def trace(x):
        return np.trace(np.array(x['stress'])[0])

    df['stress'] = df['output'].map(trace)
    plt.plot(df['factor'], df['stress'])
    plt.savefig(filename + '.pdf')


if __name__ == "__main__":
    main(sys.argv)
