from numpy import loadtxt
from pychemia.utils.periodic import atomic_symbol
from pychemia import Structure


def load(filename):
    symbols = loadtxt(filename, skiprows=2, usecols=[0], dtype='|S2', ndmin=1)
    symbols = [x.decode('utf-8') for x in symbols]

    for i in range(len(symbols)):
        if symbols[i].isdigit():
            symbols[i] = atomic_symbol(int(symbols[i]))

    positions = loadtxt(filename, skiprows=2, usecols=(1, 2, 3), ndmin=2)
    natom = len(symbols)
    periodicity = 3 * [False]
    return Structure(symbols=symbols, positions=positions, periodicity=periodicity, natom=natom)


def save(structure, filename):

    if isinstance(structure, Structure):
        sts = [structure]
    else:
        sts = structure

    wf = open(filename, 'w')
    for st in sts:
        xyz = str(st.natom) + '\n\n'
        for i in range(structure.natom):
            xyz += " %2s %15.7f %15.7f %15.7f\n" % (structure.symbols[i],
                                                    structure.positions[i, 0],
                                                    structure.positions[i, 1],
                                                    structure.positions[i, 2])
        wf.write(xyz)
    wf.close()
