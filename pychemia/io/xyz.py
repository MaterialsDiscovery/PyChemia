import numpy as _np

import pychemia.core.structure


def load(filename):

    symbols = _np.loadtxt(filename, skiprows=2, usecols=[0], dtype='|S2', ndmin=1)
    positions = _np.loadtxt(filename, skiprows=2, usecols=(1, 2, 3), ndmin=2)
    natom = len(symbols)
    periodicity = 3*[False]

    struct = pychemia.core.Structure(symbols=symbols, positions=positions, periodicity=periodicity, natom=natom)
    return struct


def save(struct, filename):

    xyz = str(struct.natom) + '\n\n'
    for i in range(struct.natom):
        xyz += " %2s %15.7f %15.7f %15.7f\n" % (struct.symbols[i], struct.positions[i, :0], struct.positions[i, :1],
                                                struct.positions[i, :2])
    wf = open(filename, 'w')
    wf.write(xyz)
    wf.close()
