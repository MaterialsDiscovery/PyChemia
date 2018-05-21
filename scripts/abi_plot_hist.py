import math
import sys

import numpy as np
from scipy.io import netcdf_file as _netcdf_file

if 'matplotlib' not in sys.modules:
    import matplotlib

    matplotlib.use('agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from pychemia.code.abinit import AbinitInput, AbiFiles
from pychemia.utils.periodic import covalent_radius
from pychemia.utils.constants import bohr_angstrom, angstrom_bohr


# 2D Plots for ABIPYTHON (Requires MATPLOTLIB)
def plot_history_energy(etotal, ekin, fcart, labels, bonds, filep):
    """
    Create a window with 2 plots, the total energy vs iterations
    and the magnitute of forces vs iterations
    """

    pp = PdfPages(filep)
    plt.ioff()

    # Page No 1
    plt.figure(figsize=(16, 10), dpi=100)

    left = 0.1
    width = 0.8
    bottom = 0.1

    if ekin is None:
        height = 0.40
    else:
        height = 0.3

    ax1 = plt.axes([left, bottom, width, height])
    ax1.plot(etotal, lw=2)
    ax1.grid(True)
    plt.ylabel("Total Energy")
    plt.xlabel("Iterations")

    if ekin is not None:
        bottom += height
        ax2 = plt.axes([left, bottom, width, height], sharex=ax1)
        ax2.semilogy(ekin, lw=2)
        ax2.grid(True)
        plt.ylabel("Ionic Kinetic Energy")
        # xticklabels = ax2.get_xticklabels()

    bottom += height
    ax3 = plt.axes([left, bottom, width, height], sharex=ax1)
    forces = fcart.transpose((1, 0, 2))
    # xticklabels = xticklabels + ax3.get_xticklabels()

    for i in range(len(forces)):
        list2 = []
        for j in range(len(forces[0])):
            list2.append(math.sqrt(sum(forces[i, j, :] * forces[i, j, :])))
        ax3.semilogy(list2, label=labels[i], lw=2)

    ax3.grid(True)
    ax3.axhline(0, color='black', lw=2)
    plt.ylabel("Forces")

    pp.savefig()

    # Page No 2
    plt.figure(figsize=(10, 16), dpi=100)

    left = 0.1
    width = 0.8
    height = 0.9 / float(len(forces))
    bottom = 0.95 - height
    plt.suptitle("Forces")

    for i in range(len(forces)):
        mag_force = []
        for j in range(len(forces[0])):
            mag_force.append(math.sqrt(sum(forces[i, j, :] * forces[i, j, :])))
        ax1 = plt.axes([left, bottom, width, height])
        ax1.semilogy(mag_force, label=labels[i], lw=2)
        ax1.grid(True)
        bottom -= height
        plt.ylabel(labels[i])
        (xmin, xmax) = ax1.get_xlim()
        if i < len(forces) - 1:
            ax1.set_xticklabels([])
        if xmax < 20:
            plt.xticks(range(int(xmax)))

    pp.savefig()

    # Page No 3
    plt.figure(figsize=(10, 16), dpi=100)

    bonds_dict = {}
    # t is time
    for t in range(len(bonds)):
        # print "Plotting bonds for time:",t
        # print "Number of bonds",len(bonds[t])
        for b in range(len(bonds[t])):

            bond = bonds[t][b]
            key = labels[bond[0]] + ' - ' + labels[bond[1]]
            if key in bonds_dict:
                bonds_dict[key].append([t, bond[2]])
            else:
                bonds_dict[key] = [[t, bond[2]]]

    left = 0.1
    width = 0.8
    nbonds = len(list(bonds_dict.keys()))
    height = 0.9 / float(nbonds)
    bottom = 0.95 - height
    plt.suptitle("Bonds")

    for i in bonds_dict:
        # print i

        ax1 = plt.axes([left, bottom, width, height])
        array = np.array(bonds_dict[i])
        ax1.plot(array[:, 0], bohr_angstrom * (array[:, 1] - array[0, 1]),
                 label='Initial=' + str(bohr_angstrom * array[0, 1]))
        ax1.grid(True)
        plt.ylabel(i)
        (xmin, xmax) = ax1.get_xlim()
        if i != list(bonds_dict.keys())[-1]:
            ax1.set_xticklabels([])
        if xmax < 20:
            plt.xticks(range(int(xmax)))
        bottom -= height
        plt.legend(loc=1)

    pp.savefig()

    pp.close()


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
    if isinstance(typat, int):
        lsttypat = [typat]
    else:
        lsttypat = typat
    if isinstance(znucl, (int, float)):
        lstznucl = [znucl]
    else:
        lstznucl = znucl

    covrad = [angstrom_bohr * covalent_radius(lstznucl[i - 1]) for i in lsttypat]
    bonds = []
    for t in range(len(xcart)):
        bonds1 = []
        for i in range(len(xcart[t])):
            for j in range(i + 1, len(xcart[t])):
                # print xcart[t][j]
                # print xcart[t][i]
                # Compute bond length between atoms i and j
                bl = math.sqrt(sum((xcart[t][j] - xcart[t][i]) ** 2))
                # print bl
                # print i,j
                # print covrad[i],covrad[j]
                if 1.35 * (covrad[i] + covrad[j]) > bl:
                    bonds1.append([i, j, bl, (xcart[t][j] - xcart[t][i]) / bl])
        bonds.append(bonds1)
    return bonds


def compute_angles(bonds, natom):
    print('natom', natom)
    itime = 0

    # Compute the number of bonds per atom
    nbonds = np.zeros(natom)
    for i in bonds[itime]:
        nbonds[i[0]] += 1
        nbonds[i[1]] += 1

    # Compute the angle for atoms with more than 1 bond
    for iatom in range(natom):
        if nbonds[iatom] > 1:
            index = []
            j = 0
            for b in bonds[itime]:
                if b[0] == iatom or b[1] == iatom:
                    index.append(j)
                j += 1
            print(index)

            for i in range(len(index)):
                for j in range(i + 1, len(index)):
                    print('i=', i)
                    print('j=', j)
                    # print nbonds


def get_history(abinitfile, dataset=""):
    if dataset == "":
        filep = abinitfile.basedir + "/" + abinitfile.files['tmpout'] + "_HIST"
    else:
        filep = abinitfile.basedir + "/" + abinitfile.files['tmpout'] + "_DS" + str(dataset) + "_HIST"
        # print filename
    # inp = AbinitInput(abinitfile.get_input_filename())

    ret = _netcdf_file(filep, 'r', mmap=False)

    return ret


def plot_history(abinitfile, dataset=""):
    history = get_history(abinitfile, dataset)
    av = AbinitInput(abinitfile.get_input_filename())

    if dataset == "":
        filep = abinitfile.basedir + "/" + abinitfile.files['tmpout'] + ".pdf"
    else:
        filep = abinitfile.basedir + "/" + abinitfile.files['tmpout'] + "_DS" + str(dataset) + ".pdf"

    xcart = history.variables['xcart'][:]
    fcart = history.variables['fcart'][:]
    # rprimd = history.variables['rprimd'][:]
    etotal = history.variables['etotal'][:]
    if 'ekin' in history.variables:
        ekin = history.variables['ekin'][:]
        if max(history.variables['ekin'][:]) == 0.0:
            ekin = None
    else:
        ekin = None

    # Getting Labels of atoms
    labels = [av.atom_name(i) for i in range(av.get_value('natom', dataset))]
    znucl = av.get_value('znucl', idtset=dataset)
    typat = av.get_value('typat', idtset=dataset)

    bonds = compute_bonds(typat, xcart, znucl)
    # natom = len(xcart[0])
    # ComputeAngles(bonds,natom)

    plot_history_energy(etotal, ekin, fcart, labels, bonds, filep)


if __name__ == '__main__':

    # Script starts from here
    if len(sys.argv) < 2:
        print('No option specified.')
        sys.exit()

    filename = sys.argv[1]
    if len(sys.argv) == 3:
        idtset = int(sys.argv[2])
    else:
        idtset = ""

    abifile = AbiFiles(filename)
    plot_history(abifile, idtset)
