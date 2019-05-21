#!/usr/bin/env python

import os
import sys
import pychemia
import argparse

if __name__ == "__main__":

    description = """Create a surface for a given structure along the given miller indices and for a given
    number of layers"""

    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-m', '--miller', nargs=3,
                        default=(0, 0, 1), metavar='N', type=int,
                        help='Miller indices (default: 0 0 1)')
    parser.add_argument('-l', '--layers',
                        default=1, metavar='N', type=int,
                        help='Number of layers (default: 1)')
    parser.add_argument('-i', '--input',
                        default='POSCAR', metavar='filename', type=str,
                        help='Input file (default: POSCAR)')
    parser.add_argument('-o', '--output',
                        default='POSCARsurf', metavar='filename', type=str,
                        help='Output file (default: POSCARsurf)')

    args = parser.parse_args()

    inputfile = args.input
    outputfile = args.output
    miller = tuple(args.miller)
    layers = args.layers

    if not os.path.isfile(inputfile):
        print('File not found: ', input)
        sys.exit(1)

    print('INPUT:  ', inputfile)
    print('OUTPUT: ', outputfile)
    print('MILLER: ', miller)
    print('LAYERS: ', layers)

    hh = int(miller[0])
    kk = int(miller[1])
    ll = int(miller[2])

    structure = pychemia.code.vasp.read_poscar(inputfile)
    surf = pychemia.analysis.rotate_along_indices(structure, hh, kk, ll, layers)

    print(surf)

    pychemia.code.vasp.write_poscar(surf, outputfile)

    symmetry = pychemia.crystal.CrystalSymmetry(structure)
    print('SpaceGroup (Original): ', symmetry.number())

    symmetry = pychemia.crystal.CrystalSymmetry(surf)
    print('SpaceGroup (Final)   : ', symmetry.number())
