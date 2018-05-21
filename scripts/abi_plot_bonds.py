import math
import sys

import numpy as np

if 'matplotlib' not in sys.modules:
    import matplotlib

    matplotlib.use('agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from pychemia.code.abinit import AbinitInput, AbiFiles
from pychemia.utils.periodic import covalent_radius
from pychemia.utils.constants import bohr_angstrom


def compute_bonds(typat, xcart, znucl):
    """
    Compute the bond lengths of all the atoms
    inside the unitary box (NEEDS EXTENSION TO
    OUTSIDE THE BOX)

    :param typat: (int, list) Type of atoms
    :param xcart: (numpy.ndarray) Cartesian positions
    :param znucl: (int, list) Atomic number for atoms in typat
    :return:
    """
    print(typat)
    print(xcart)
    print(znucl)

    npxcart = np.array(xcart).reshape((-1, 3))
    if isinstance(typat, int):
        lsttypat = [typat]
    else:
        lsttypat = typat
    if isinstance(znucl, (int, float)):
        lstznucl = [znucl]
    else:
        lstznucl = znucl
    covrad = [covalent_radius(lstznucl[iatom - 1]) for iatom in lsttypat]

    bonds = []
    for iatom in range(len(npxcart)):
        for jatom in range(iatom + 1, len(npxcart)):
            # Compute bond length between atoms i and j
            bl = math.sqrt(sum((npxcart[jatom] - npxcart[iatom]) ** 2))
            if 1.35 * (covrad[iatom] + covrad[jatom]) > bl:
                bonds.append([iatom, jatom, bl, (npxcart[jatom] - npxcart[iatom]) / bl])
            else:
                print("small distance: %s" % bl)

    return bonds


def save_bonds(allbonds, abivar):
    wf = open("bonds.txt", 'w')

    for ii in range(len(allbonds[0])):  # Number of bonds
        label = abivar.atom_name(allbonds[0][ii][0]) + ' ' + abivar.atom_name(allbonds[0][ii][1])
        wf.write(label)
        for jj in range(len(allbonds)):  # Number of sets
            wf.write(' ' + str(round(bohr_angstrom * allbonds[jj][ii][2], 4)).ljust(10))  # Bond length
        wf.write('\n')
    wf.close()


def get_all_bonds(abinitfile, dtset=''):
    allbonds = []
    index = 0
    for abf in abinitfile:
        if dtset == '':
            idt = ''
        else:
            idt = int(dtset[index])
            index += 1
        filep = abf.basedir + "/" + abf.files['tmpout'] + "_OUT.nc"

        out = AbinitInput(filep)
        
        print(out.variables)

        xcart = out.get_value('xcart', idtset=idt)
        znucl = out.get_value('znucl', idtset=idt)
        typat = out.get_value('typat', idtset=idt)

        bonds = compute_bonds(typat, xcart, znucl)
        allbonds.append(bonds)

    abivar = AbinitInput(abinitfile[0].get_input_filename())
    save_bonds(allbonds, abivar)
    return allbonds


def plot_bonds(listabifile, listidtset):
#    print(listabifile)
#    print(listidtset)
    pp = PdfPages("bonds.pdf")
    plt.ioff()

    if not listidtset:
        allbonds = get_all_bonds(listabifile)
    else:
        allbonds = get_all_bonds(listabifile, listidtset)
    print(allbonds)

    plt.figure(figsize=(32, 20), dpi=100)
    plt.subplots_adjust(left=0.05, bottom=0.05, right=0.95, top=0.95,
                        wspace=None, hspace=None)

    abivar = AbinitInput(listabifile[0].get_input_filename())

    nbond = len(allbonds[0])
    n = 0
    for ibond in allbonds:
        if len(ibond) != nbond:
            print('ERROR: Number of bonds changed')
            if not listidtset:
                print(listabifile[n].basedir, 'contains', len(ibond), 'bonds')
            else:
                print(listabifile[n].filename + ':' + listidtset[n], 'contains', len(ibond), 'bonds')
                n += 1

    iplot = 1
    for ibond in range(nbond):
        plt.subplot(5, 4, iplot)
        y = [bohr_angstrom * seti[ibond][2] for seti in allbonds]
        print(abivar.atom_name(seti[ibond][0]))
        print(abivar.atom_name(seti[ibond][1]))
        label = abivar.atom_name(seti[ibond][0]) + ' ' + abivar.atom_name(seti[ibond][1])
        plt.plot(y, label=label)
        plt.plot(y, 'ro')
        iplot += 1
        plt.legend()

    pp.savefig()
    plt.clf()

    for ibond in range(nbond):
        y = [bohr_angstrom * seti[ibond][2] for seti in allbonds]
        label = abivar.atom_name(seti[ibond][0]) + ' ' + abivar.atom_name(seti[ibond][1])
        plt.plot(y, label=label)
        plt.plot(y, 'ro')

        plt.text(0.09, y[0] + 0.001, label, size='small')
        iplot += 1
        plt.legend()

    pp.savefig()
    pp.close()


if __name__ == '__main__':

    narg = len(sys.argv) - 1

    list_abifile = []
    list_idtset = []
    dtsetall = False
    for i in range(narg):
        filename = sys.argv[i + 1]
        if ':' in filename:
            idtset = filename.split(':')[1]
            if idtset == 'all':
                dtsetall = True
            else:
                list_idtset.append(idtset)
            filename = filename.split(':')[0]

        abifile = AbiFiles(filename)
        if dtsetall:
            av = AbinitInput(abifile.get_input_filename())
            keys = av.get_dtsets_keys()
#            print(keys)
            for j in keys:
                list_abifile.append(abifile)
                list_idtset.append(str(j))

        else:
            list_abifile.append(abifile)

    plot_bonds(list_abifile, list_idtset)
