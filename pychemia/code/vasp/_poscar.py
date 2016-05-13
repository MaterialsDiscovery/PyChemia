"""
Routines to read and write POSCAR file
"""

import os
import numpy as _np
import pychemia


def read_poscar(path='POSCAR'):
    """
    Load a POSCAR file and return a pychemia structure object

    :param path: (str) Filename of the POSCAR to read
    :return:
    """

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
        print("POSCAR path not found")
        return

    # Reading the POSCAR file
    rf = open(poscarfile, 'r')
    comment = rf.readline().strip()
    latconst = float(rf.readline())
    newcell = _np.zeros((3, 3))

    newcell[0, :] = latconst * _np.array([float(x) for x in rf.readline().split()])
    newcell[1, :] = latconst * _np.array([float(x) for x in rf.readline().split()])
    newcell[2, :] = latconst * _np.array([float(x) for x in rf.readline().split()])

    line = rf.readline()
    species = None

    try:
        natom_per_species = _np.array([int(x) for x in line.split()])
        # print 'Old Format'
    except ValueError:
        # print 'New format'
        species = [x for x in line.split()]
        line = rf.readline()
        natom_per_species = _np.array([int(x) for x in line.split()])

    natom = _np.sum(natom_per_species)

    if species is None:
        if os.path.isfile(potcarfile):
            species = get_species(potcarfile)
        else:
            print(""" ERROR: The POSCAR does not contain information about the species present on the structure
            You can set a consistent POTCAR along the POSCAR or
            modify your POSCAR by adding the atomic symbol on the sixth line of the file""")
            return None

    if not species:
        print('No information about species')
        raise ValueError()

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

    pos = []
    for i in range(natom):
        pos += [float(x) for x in rf.readline().split()[:3]]
    pos = _np.array(pos).reshape((-1, 3))

    if kmode == 'Cartesian':
        return pychemia.Structure(cell=newcell, symbols=symbols, reduced=pos, comment=comment)
    else:
        return pychemia.Structure(cell=newcell, symbols=symbols, reduced=pos, comment=comment)


def write_poscar(structure, filepath='POSCAR', newformat=True):
    """
    Takes an structure from pychemia and save the file
    POSCAR for VASP.

    :param structure: (pychemia.Structure) Structure to write POSCAR
    :param filepath: (str) Filename of POSCAR file to create
    :param newformat: (bool) If the new VASP format is used to create the POSCAR
    """
    ret = ''
    comp = structure.get_composition()
    for i in comp.species:
        ret += ' ' + i
    ret += '\n'
    ret += '1.0\n'
    for i in range(3):
        ret += ' %20.16f %20.16f %20.16f\n' % tuple(structure.cell[i])
    if newformat:
        for i in comp.species:
            ret += ' ' + i
        ret += '\n'
    for i in comp.values:
        ret += ' ' + str(i)
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
        if 'PAW' in line and 'PAW' == line.split()[0].strip() and 'radial' not in line:
            species.append(line.split()[1].split('_')[0])

    return species


def write_potcar(structure, filepath='POTCAR', pspdir='potpaw_PBE', options=None, pspfiles=None):
    comp = structure.get_composition()
    ret = ''
    psppath = os.getenv('HOME') + '/.vasp/PP-VASP/' + pspdir
    if not os.path.exists(psppath):
        raise ValueError("The path for VASP Pseudo-potentials does not exists: " + psppath)

    if pspfiles is None:
        pspfiles = []
        for i in comp.species:
            if options is not None and i in options:
                if isinstance(options[i], basestring):
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
                    else:
                        pass
                        # print pspfile, 'is not present...'
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
