"""
Routines to read and write POSCAR file
"""

import os

import numpy as _np

import pychemia


def load_POSCAR(path):
    """
    Load a POSCAR file and return a pychemia crystal object
    """

    if os.path.isfile(path):
        filename = path
    elif os.path.isdir(path) and os.path.isfile(path + '/POSCAR'):
        filename = path + '/POSCAR'
    else:
        print("POSCAR path not found")
        return

    crystal = pychemia.geometry.Structure()

    # Reading the POSCAR file
    rf = open(filename, 'r')
    crystal.comment = rf.readline().strip()
    latconst = float(rf.readline())
    crystal.cell[0, :] = _np.array([float(x) for x in rf.readline().split()])
    crystal.cell[1, :] = _np.array([float(x) for x in rf.readline().split()])
    crystal.cell[2, :] = _np.array([float(x) for x in rf.readline().split()])

    # The call to add_atom naturally will increase the
    # internal variable crystal.natom
    natom_per_species = _np.array([int(x) for x in rf.readline().split()])
    natom = _np.sum(natom_per_species)

    species = get_species(path + '/POTCAR')

    symbols = []
    for i in range(len(natom_per_species)):
        numspe = natom_per_species[i]
        for j in range(numspe):
            symbols.append(species[i])

    mode = rf.readline()
    if mode[0].lower() in ['c', 'k']:
        kmode = 'Cartesian'
    else:
        kmode = 'Direct'

    for i in range(natom):
        pos = [float(x) for x in rf.readline().split()]
        crystal.add_atom(symbols[i], pos)

    return crystal


def save_POSCAR(crystal, path):
    """
    Takes an object crystal from pychemia and save the file
    POSCAR for VASP.
    """
    pass


def get_species(path):
    species = []
    rf = open(path, 'r')
    for line in rf.readlines():
        if 'PAW_PBE' in line and 'PAW_PBE' == line.split()[0].strip():
            species.append(line.split()[1].split('_')[0])

    return species
