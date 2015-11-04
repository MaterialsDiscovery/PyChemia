import numpy as _np
import pychemia.core.structure


def load(filename):
    symbols = _np.loadtxt(filename, skiprows=2, usecols=[0], dtype='|S2', ndmin=1)
    positions = _np.loadtxt(filename, skiprows=2, usecols=(1, 2, 3), ndmin=2)
    natom = len(symbols)
    periodicity = 3 * [False]

    return pychemia.core.Structure(symbols=symbols, positions=positions, periodicity=periodicity, natom=natom)


def save(structure, filename):
    xyz = str(structure.natom) + '\n\n'
    for i in range(structure.natom):
        xyz += " %2s %15.7f %15.7f %15.7f\n" % (structure.symbols[i],
                                                structure.positions[i, 0],
                                                structure.positions[i, 1],
                                                structure.positions[i, 2])
    wf = open(filename, 'w')
    wf.write(xyz)
    wf.close()
