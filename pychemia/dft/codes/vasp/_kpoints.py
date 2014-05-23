"""
Routines to read and write KPOINTS file
"""

__author__ = "Guillermo Avendano-Franco"
__copyright__ = "Copyright 2014"
__version__ = "0.1"
__maintainer__ = "Guillermo Avendano-Franco"
__email__ = "gtux.gaf@gmail.com"
__status__ = "Development"
__date__ = "March 16, 2014"

import os

import numpy as _np

import pychemia


def load_KPOINTS(path):
    """
    Load the file KPOINTS in the directory 'path' or
    read directly the file 'path' and return a kpoints
    object for pychemia
    """
    if os.path.isfile(path):
        filename = path
    elif os.path.isdir(path) and os.path.isfile(path + '/KPOINTS'):
        filename = path + '/KPOINTS'
    else:
        print("KPOINTS path not found")
        return

    kp = pychemia.dft.KPoints(nkpt=None)

    # Reading the KPOINTS file
    rf = open(filename, 'r')
    kp.comment = rf.readline()
    print(kp.comment)
    nkpt = int(rf.readline())
    mode = rf.readline()
    if mode[0].lower() in ['c', 'k']:
        kmode = 'Cartesian'
    else:
        kmode = 'Reciprocal'

    line = _np.array([float(x) for x in rf.readline().split()])
    pos = line[:3]
    wgt = line[3]
    kp.add_kpt(pos, wgt)

    return kp


def save_KPOINTS(kp, path):
    """
    Takes an object kpoints from pychemia and
    save the file KPOINTS in the directory 'path' or
    save the file 'path' as a VASP KPOINTS file
    """

    if os.path.isdir(path):
        filename = path + '/KPOINTS'
    else:
        filename = path

    wf = open(filename, 'w')

    wf.write(kp.comment + '\n')
    wf.write(kp.nkpt + '\n')
    for i in range(kp.nkpt):
        wf.write(" %15.7f %15.7f %15.7f %20.7f\n"
                 % (kp.kpts[i, 0],
                    kp.kpts[i, 1],
                    kp.kpts[i, 2],
                    kp.wgts[i]))

    wf.close()

