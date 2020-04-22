# coding: utf-8
# Copyright (c) Guillermo Avendano-Franco
# Distributed under the terms of the MIT License.

import os
from re import findall
from numpy import zeros, array, sum
from pychemia import Structure
from pychemia.utils.periodic import atomic_symbols
from itertools import groupby

"""
Routines to read and write POSCAR file
"""

__author__ = "Guillermo Avendano-Franco"
__copyright__ = "Copyright 2017, PyChemia Project"
__version__ = "0.2"
__maintainer__ = "Guillermo Avendano-Franco"
__email__ = "gufranco@mail.wvu.edu"
__date__ = "Feb 20, 2017"


def read_poscar(path='POSCAR'):
    """
    Load a POSCAR file and return a pychemia structure object

    :param path: (str) Path to a POSCAR file or a directory where a file named 'POSCAR' is located
    :return:
    """
    # The argument 'path' could refer to a POSCAR file or the directory where the file 'POSCAR' exists
    if os.path.isfile(path):
        poscarfile = path
        if os.path.dirname(path) != '':
            potcarfile = os.path.dirname(path) + os.sep + 'POTCAR'
        else:
            potcarfile = 'POTCAR'
    elif os.path.isdir(path) and os.path.isfile(path + os.sep + 'POSCAR'):
        poscarfile = path + os.sep + 'POSCAR'
        potcarfile = path + os.sep + 'POTCAR'
    else:
        raise ValueError("ERROR: No POSCAR file found on %s" % path)

    # Reading the POSCAR file
    rf = open(poscarfile, 'r')
    comment = rf.readline().strip()
    latconst = float(rf.readline().split()[0])
    newcell = zeros((3, 3))

    newcell[0, :] = latconst * array([float(x) for x in rf.readline().split()[:3]])
    newcell[1, :] = latconst * array([float(x) for x in rf.readline().split()[:3]])
    newcell[2, :] = latconst * array([float(x) for x in rf.readline().split()[:3]])

    line = rf.readline()
    species = None

    # This is the old format, the only way of knowing which species are refering is by
    # reading the POTCAR file
    natom_per_species = array([int(x) for x in line.split() if x.isdigit()])

    # Check if the file is the new format, in such case this line contains the
    # atomic symbols of the atoms and the next line the number of atoms of each
    # species. The new format makes a POSCAR self-contained to create a structure
    if len(natom_per_species) == 0:
        species = [x for x in line.split() if x in atomic_symbols]
        line = rf.readline()
        natom_per_species = array([int(x) for x in line.split() if x.isdigit()])

    natom = sum(natom_per_species)

    if species is None:
        comment_species = [x for x in comment.split() if x in atomic_symbols]
        if os.path.isfile(potcarfile):
            species = get_species(potcarfile)
        elif len(comment_species) == len(natom_per_species):
            species = comment_species
        else:
            raise ValueError("ERROR: The POSCAR does not have information about the species present on the structure\n"
                             + "You can set a consistent POTCAR along the POSCAR or modify your POSCAR to the\n"
                             + "new format by adding the atomic symbol(s) on the sixth line of the file")

    symbols = []
    for i in range(len(natom_per_species)):
        for j in range(natom_per_species[i]):
            symbols.append(species[i])

    mode = rf.readline()
    if mode[0].lower() in ['c', 'k']:
        kmode = 'Cartesian'
    else:
        kmode = 'Direct'

    pos = []
    for i in range(natom):
        pos += [float(x) for x in rf.readline().split()[:3]]
    pos = array(pos).reshape((-1, 3))

    if kmode == 'Cartesian':
        return Structure(cell=newcell, symbols=symbols, positions=pos, comment=comment)
    else:
        return Structure(cell=newcell, symbols=symbols, reduced=pos, comment=comment)


def write_poscar(structure, filepath='POSCAR', newformat=True, direct=True, comment=None, heterostructure=False):
    """
    Takes an structure from pychemia and save the file
    POSCAR for VASP.

    :param comment: Optional comment to the first line of the POSCAR
    :param structure: (pychemia.Structure) Structure to write POSCAR
    :param filepath: (str) Filename of POSCAR file to create
    :param newformat: (bool) If the new VASP format is used to create the POSCAR
    :param direct: (bool) If True, use reduced coordinates. If False, use cartesian coordinates (default: True)
    """
    comp = structure.get_composition()

    # If heterostructure is true it will keep the repeating order found
    # in the POSCAR. 
    # Added by Uthpala on Apr 20th, 2020.
    if heterostructure:
        species = [i[0] for i in groupby(structure.symbols)]
    else:
        species = get_species_list(structure)
    species_count = [len(list(group)) for key, group in groupby(structure.symbols)]   

    ret = ''
    if comment is None:
        for i in species:
            ret += ' ' + i
    else:
        ret += comment.strip()
    ret += '\n'
    ret += '1.0\n'
    for i in range(3):
        ret += ' %20.16f %20.16f %20.16f\n' % tuple(structure.cell[i])

    if newformat:
        for i in species:
            ret += ' ' + i
        ret += '\n'
    for icount, i in enumerate(species):
        if heterostructure:
           ret += ' ' + str(species_count[icount]) 
        else:
            ret += ' ' + str(comp.composition[i])
    ret += '\n'

    if direct:
        ret += 'Direct\n'
        for i in range(structure.natom):
            ret += ' %20.16f %20.16f %20.16f\n' % tuple(structure.reduced[i])
    else:
        ret += 'Cartesian\n'
        for i in range(structure.natom):
            ret += ' %20.16f %20.16f %20.16f\n' % tuple(structure.positions[i])

    wf = open(filepath, 'w')
    wf.write(ret)
    wf.close()


def get_species_list(structure):
    while True:
        species = []
        for i in structure.symbols:
            if i not in species:
                species.append(i)
        if len(species) == len(structure.species):
            break
        else:
            structure.sort_sites()
    return species


def get_species(path):
    species = []
    rf = open(path, 'r')
    for line in rf.readlines():
        if 'PAW_PBE' in line and 'PAW_PBE' == line.split()[0].strip():
            species.append(line.split()[1].split('_')[0])
        if 'PAW' in line and 'PAW' == line.split()[0].strip() and 'radial' not in line:
            species.append(line.split()[1].split('_')[0])

    return species


def write_potcar(structure, filepath='POTCAR', pspdir='potpaw_PBE', options=None, pspfiles=None, basepsp=None, heterostructure=False):

    # If heterostructure is true it will keep the repeating order found
    # in the POSCAR. 
    # Added by Uthpala on Apr 20th, 2020.
    if heterostructure:
        species = [i[0] for i in groupby(structure.symbols)]
    else:
        species = get_species_list(structure)

    ret = ''
    if basepsp is not None:
        psppath = os.path.abspath(basepsp) + os.sep + pspdir
    else:
        psppath = os.getenv('HOME') + '/.vasp/PP-VASP/' + pspdir
    if not os.path.exists(psppath):
        raise ValueError("The path for VASP Pseudo-potentials does not exists: " + psppath)

    if pspfiles is None:
        pspfiles = []
        for i in species:
            if options is not None and i in options:
                if isinstance(options[i], str):
                    pspfile = psppath + os.sep + i + '_' + options[i] + '/POTCAR'
                elif isinstance(options[i], list):
                    for j in options[i]:
                        pspfile = psppath + os.sep + i + '_' + options[i][j] + '/POTCAR'
                        if os.path.isfile(psppath):
                            break

            else:
                for j in ['', '_sv']:
                    pspfile = psppath + os.sep + i + j + '/POTCAR'
                    if os.path.isfile(pspfile):
                        break
            if not os.path.isfile(pspfile):
                raise ValueError("File not found : " + pspfile)
            pspfiles.append(pspfile)

    for pspfile in pspfiles:
        rf = open(pspfile)
        ret += rf.read()
        rf.close()
    wf = open(filepath, 'w')
    wf.write(ret)
    wf.close()
    return pspfiles


def get_potcar_info(filename='POTCAR'):
    rf = open(filename)
    data = rf.read()
    pairs = findall('([\w ]*)=([.\d ]*)', data)
    ret = {}
    for i in pairs:
        print(i)
        if i[1].strip() == '':
            pass
        elif i[1].strip().isdigit():
            ret[i[0].strip()] = int(i[1])
        elif len(i[1].split()) > 1:
            ret[i[0].strip()] = [float(x) for x in i[1].split()]
        else:
            ret[i[0].strip()] = float(i[1])

    return ret
