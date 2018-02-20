#!/usr/bin/env python

import logging
import os
import sys
import json
import getopt
import pychemia
from pychemia.code.vasp.task import ConvergenceKPointGrid, IonRelaxation, ConvergenceCutOffEnergy
from pychemia.utils.computing import get_float, get_int

logging.basicConfig(level=logging.DEBUG)

try:
    from pychemia.symm import symmetrize
except ImportError:
    def symmetrize(x):
        return x


def cleaner():
    files = ['DOSCAR', 'EIGENVAL', 'IBZKPT', 'OSZICAR', 'PCDAT', 'PROCAR', 'WAVECAR', 'XDATCAR']
    [os.remove(ifile) for ifile in files if os.path.exists(ifile)]


def usage(name):
    print("""
NAME
    %s

DESCRIPTION
    Vasp relaxator with convergence in K-Point Grid and Cut-Off Energy (ENCUT)

OPTIONS

    --help, -h
        Return information on the options and use of this script

    --structure, -s <string> (Default: 'POSCAR')
        Structure for which the strength will be computed

    --output, -o <string> (Default: 'pcm_results.json')
        File in JSON format where the results are stored

    --nparal, -n <int> (Default: 2)
        Number of MPI parallel processes for VASP

    --binary, -b <string> (Default: 'vasp')
        Path to the VASP executable

    --energy_tol, -e <float> (Default: 1E-3)
        Energy tolerance for three succesive values of the Kpoint grid or CutOff energy

    --target_forces, -t <float> (Default: 1E-3)
        Target forces for internal relaxation
""" % os.path.basename(name))


def main(argv):
    try:
        opts, args = getopt.getopt(argv[1:], "hs:o:n:b:e:t:", ["help", "structure=", "output=", "nparal=",
                                                               "binary=", "energy_tol=", "target_forces="])
    except getopt.GetoptError:
        usage(argv[0])
        sys.exit(2)

    if len(opts) == 0:
        usage(argv[0])
        sys.exit(2)

    # Default Values
    workdir = '.'
    structure_file = 'POSCAR'
    output_file = 'pcm_results.json'
    nparal = 2
    energy_tol = 1E-3
    target_forces = 1E-3
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
        elif opt in ("-t", "--target_forces"):
            target_forces = get_float(arg)
        elif opt in ("-s", "--structure"):
            structure_file = arg

    structure = pychemia.structure_from_file(structure_file)
    if structure is None:
        print(" ERROR: Invalid structure, no structure could be obtained from '%s'" % structure_file)
        sys.exit(2)

    if structure_file == 'POSCAR':
        os.rename('POSCAR', 'POSCAR_original')

    print("\n PyChemia VASP Relaxator")
    print(" =======================\n")
    print(" VASP binary         : ", binary)
    print(" Energy tolerance    : ", energy_tol)
    print(" Target forces       : ", target_forces)
    print(" MPI number of procs : ", nparal)
    print(" Structure           :\n")
    print(structure)

    wf = open(output_file, 'w')
    data = {'input': {'binary': binary,
                      'energy_tol': energy_tol,
                      'target_forces': target_forces,
                      'nparal': nparal,
                      'structure': structure.to_dict}}
    json.dump(data, wf)
    wf.close()

    # First Round (Relaxing the original structure)
    print('\nFirst Round')
    print('===========')

    cleaner()
    print('\nConvergence of Cut-off Energy')
    print('-----------------------------\n')
    ce = ConvergenceCutOffEnergy(structure, energy_tolerance=energy_tol)
    if os.path.isfile('convergence_encut.json'):
        print('A previous convergence study was found, loading...')
        ce.load()
    if not ce.is_converge:
        ce.run(nparal)
        ce.save()
        ce.plot()
    encut = ce.best_encut

    data['output'] = {'1R_ENCUT': encut}
    wf = open(output_file, 'w')
    json.dump(data, wf)
    wf.close()

    cleaner()
    print('\nConvergence of K-Point Grid')
    print('---------------------------\n')
    ck = ConvergenceKPointGrid(structure, encut=encut, energy_tolerance=energy_tol)
    if os.path.isfile('convergence_kpoints.json'):
        print('A previous convergence study was found, loading...')
        ce.load()
    if not ce.is_converge:
        ce.run(nparal)
        ce.save()
        ce.plot()
    kp = ck.best_kpoints

    data['output'] = {'1R_KPOINTS': list(kp.grid)}
    print(data)
    wf = open(output_file, 'w')
    json.dump(data, wf)
    wf.close()

    os.rename('convergence_encut.json', 'convergence_encut_phase1.json')
    os.rename('convergence_encut.pdf', 'convergence_encut_phase1.pdf')
    os.rename('convergence_kpoints.json', 'convergence_kpoints_phase1.json')
    os.rename('convergence_kpoints.pdf', 'convergence_kpoints_phase1.pdf')

    cleaner()
    print('\nIonic Relaxation')
    print('----------------\n')
    vr = IonRelaxation(structure=structure, encut=encut, kp_grid=kp.grid, workdir=workdir,
                       target_forces=10 * target_forces)
    vr.run(nparal)

    structure = vr.get_final_geometry()
    structure.save_json(workdir + os.sep + 'structure_phase1.json')

    data['output'] = {'1R_structure': structure.to_dict}
    wf = open(output_file, 'w')
    json.dump(data, wf)
    wf.close()

    # Second Round (Symetrize structure and redo convergences)
    st = symmetrize(structure)

    print('\nSecond Round')
    print('============')

    cleaner()
    print('\nConvergence of K-Point Grid')
    print('---------------------------\n')
    ck = ConvergenceKPointGrid(st, encut=encut, energy_tolerance=energy_tol, recover=True)
    ck.run(nparal)
    ck.save()
    ck.plot()
    kp = ck.best_kpoints

    data['output'] = {'2R_KPOINTS': list(kp.grid)}
    wf = open(output_file, 'w')
    json.dump(data, wf)
    wf.close()

    cleaner()
    print('\nConvergence of Cut-off Energy')
    print('-----------------------------\n')
    ce = ConvergenceCutOffEnergy(st, energy_tolerance=energy_tol, kpoints=kp)
    ce.run(nparal)
    ce.save()
    ce.plot()
    encut = ce.best_encut

    data['output'] = {'2R_ENCUT': encut}
    wf = open(output_file, 'w')
    json.dump(data, wf)
    wf.close()

    os.rename('convergence_encut.json', 'convergence_encut_phase2.json')
    os.rename('convergence_encut.pdf', 'convergence_encut_phase2.pdf')
    os.rename('convergence_kpoints.json', 'convergence_kpoints_phase2.json')
    os.rename('convergence_kpoints.pdf', 'convergence_kpoints_phase2.pdf')

    cleaner()
    print('\nIonic Relaxation')
    print('----------------\n')
    vr = IonRelaxation(structure=st, workdir='.', encut=encut, kp_grid=kp.grid, target_forces=target_forces)
    vr.run(nparal)

    structure = vr.get_final_geometry()
    structure.save_json(workdir + os.sep + 'structure_phase2.json')

    data['output'] = {'2R_structure': structure.to_dict}

    wf = open(output_file, 'w')
    json.dump(data, wf)
    wf.close()

    cleaner()


if __name__ == "__main__":
    main(sys.argv)
