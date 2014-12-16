"""
Routines to read and write POSCAR file
"""

import os
import numpy as _np

import pychemia


def read_poscar(path):
    """
    Load a POSCAR file and return a pychemia structure object
    """

    if os.path.isfile(path):
        poscarfile = path
        if os.path.dirname(path) != '':
            potcarfile = os.path.dirname(path)+os.sep+'POTCAR'
        else:
            potcarfile = 'POTCAR'
    elif os.path.isdir(path) and os.path.isfile(path + os.sep + 'POSCAR'):
        poscarfile = path + os.sep + 'POSCAR'
        potcarfile = path + os.sep + 'POTCAR'
    else:
        print("POSCAR path not found")
        return

    structure = pychemia.core.Structure()

    # Reading the POSCAR file
    rf = open(poscarfile, 'r')
    structure.comment = rf.readline().strip()
    latconst = float(rf.readline())
    newcell = _np.zeros((3, 3))

    newcell[0, :] = latconst * _np.array([float(x) for x in rf.readline().split()])
    newcell[1, :] = latconst * _np.array([float(x) for x in rf.readline().split()])
    newcell[2, :] = latconst * _np.array([float(x) for x in rf.readline().split()])

    structure.set_cell(newcell)

    # The call to add_atom naturally will increase the
    # internal variable crystal.natom
    line = rf.readline()
    species = None
    try:
        natom_per_species = _np.array([int(x) for x in line.split()])
    except ValueError:
        species = [x for x in line.split()]
        line = rf.readline()
        natom_per_species = _np.array([int(x) for x in line.split()])

    natom = _np.sum(natom_per_species)

    if species is None:
        species = get_species(potcarfile)

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
        if kmode == 'Cartesian':
            structure.add_atom(symbols[i], pos, option='cartesian')
        elif kmode == 'Direct':
            structure.add_atom(symbols[i], pos, option='reduced')

    return structure


def write_poscar(structure, filepath='POSCAR'):
    """
    Takes an structure from pychemia and save the file
    POSCAR for VASP.
    """
    ret = ''
    comp = structure.get_composition()
    for i in comp.species:
        ret += ' '+i
    ret += '\n'
    ret += '1.0\n'
    for i in range(3):
        ret += ' %20.16f %20.16f %20.16f\n' % tuple(structure.cell[i])
    for i in comp.values:
        ret += ' '+str(i)
    ret += '\n'
    ret += 'Direct\n'
    for i in range(structure.natom):
        ret += ' %20.16f %20.16f %20.16f\n' % tuple(structure.reduced[i])
    wf = open(filepath, 'w')
    wf.write(ret)
    wf.close()


def get_species(path):
    species = []
    rf = open(path, 'r')
    for line in rf.readlines():
        if 'PAW_PBE' in line and 'PAW_PBE' == line.split()[0].strip():
            species.append(line.split()[1].split('_')[0])

    return species


def write_potcar(structure, filepath='POTCAR', pspdir='potpaw_PBE', options=None, pspfiles=None):

    comp = structure.get_composition()
    ret = ''
    psppath = os.getenv('HOME') + '/.vasp/PP-VASP/'+pspdir
    if not os.path.exists(psppath):
        raise ValueError("The path for VASP Pseudo-potentials does not exists: "+psppath)

    if pspfiles is None:
        pspfiles = []
        for i in comp.species:
            if options is not None and i in options:
                if isinstance(options[i], basestring):
                    pspfile = psppath+os.sep+i+'_'+options[i]+'/POTCAR'
                elif isinstance(options[i], list):
                    for j in options[i]:
                        pspfile = psppath+os.sep+i+'_'+options[i][j]+'/POTCAR'
                        if os.path.isfile(psppath):
                            break

            else:
                for j in ['', '_sv']:
                    pspfile = psppath+os.sep+i+j+'/POTCAR'
                    if os.path.isfile(pspfile):
                        break
                    else:
                        print pspfile, 'is not present...'
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