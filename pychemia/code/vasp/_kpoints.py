import os
import numpy as _np
import pychemia

__author__ = "Guillermo Avendano-Franco"
__copyright__ = "Copyright 2016"
__version__ = "0.1"
__maintainer__ = "Guillermo Avendano-Franco"
__email__ = "gtux.gaf@gmail.com"
__status__ = "Development"
__date__ = "May 16, 2016"


def read_kpoints(path='KPOINTS'):
    """
    Load the file KPOINTS in the directory 'path' or
    read directly the file 'path' and return a kpoints
    object for pychemia

    :param path: (str) File path for KPOINTS file
    :return:
    """
    if os.path.isfile(path):
        filename = path
    elif os.path.isdir(path) and os.path.isfile(path + '/KPOINTS'):
        filename = path + '/KPOINTS'
    else:
        print("KPOINTS path not found")
        return

    # Reading the KPOINTS file
    rf = open(filename, 'r')
    comment = rf.readline()
    del comment
    nkpt = int(rf.readline())
    mode = rf.readline()
    if nkpt > 0:
        if mode[0].lower() in ['c', 'k']:
            kmode = 'Cartesian'
        else:
            kmode = 'Reciprocal'
        kp = pychemia.crystal.KPoints(kmode=kmode)
        for i in range(nkpt):
            line = _np.array([float(x) for x in rf.readline().split()])
            pos = line[:3]
            wgt = line[3]
            kp.add_kpt(pos, wgt)

    else:
        if mode[0].lower() in ['g']:
            kmode = 'Gamma'
        elif mode[0].lower() in ['m']:
            kmode = 'Monkhorst-pack'
        else:
            raise ValueError("Kpoints mode must be 'Gamma' or 'Monkhorst-pack'")
        kp = pychemia.crystal.KPoints(kmode=kmode)
        line = _np.array([int(x) for x in rf.readline().split()])
        grid = line[:3]
        try:
            line = _np.array([float(x) for x in rf.readline().split()])
            shift = line[:3]
        except ValueError:
            shift = _np.zeros(3)

        kp.set_grid(grid, shift)

    return kp


def write_kpoints(kp, filepath='KPOINTS'):
    """
    Takes an object kpoints from pychemia and
    save the file KPOINTS in the directory 'path' or
    save the file 'path' as a VASP KPOINTS file

    :param kp: Kpoints object
    :param filepath: (str) Filename where the KPOINTS file is created
    """

    if os.path.isdir(filepath):
        filename = filepath + '/KPOINTS'
    else:
        filename = filepath

    wf = open(filename, 'w')

    wf.write('Automatic mesh\n')
    if kp.kmode == 'cartesian' or kp.kmode == 'reciprocal':
        wf.write(str(kp.nkpt) + '\n')
        wf.write(kp.kmode.title() + '\n')
        for i in range(kp.nkpt):
            wf.write(" %15.7f %15.7f %15.7f %20.7f\n"
                     % (kp.kpoints_list[i, 0],
                        kp.kpoints_list[i, 1],
                        kp.kpoints_list[i, 2],
                        kp.weights[i]))
    elif kp.kmode == 'gamma' or kp.kmode == 'monkhorst-pack':
        wf.write('0\n')
        wf.write(kp.kmode.title() + '\n')
        wf.write(" %7d %7d %7d\n" % tuple(kp.grid))
        wf.write(" %7.4f %7.4f %7.4f\n" % tuple(kp.shifts))

    wf.close()
