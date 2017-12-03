"""
Useful routines to get physical data from octopus TD runs
"""

import os
from math import acos
import scipy.io.netcdf
from numpy import linalg, zeros, loadtxt, apply_along_axis, newaxis, polyfit, poly1d, dot, array
import pychemia


def get_last_iteration(dirname):
    """
    Get the last "td.xxxxxxx" found in dirname

    :param dirname: (str) Directory to take last information
    :return:
    """
    return sorted([int(x[3:]) for x in os.listdir(dirname) if x[:3] == 'td.' and x[-1].isdigit()])[-1]


def hirshfeld(dirname, iteration=None, spin=None):
    """
    Get the Hirshfeld charge analysis for all the atoms
    comparing for a given iteration and spin.

    Args:
        dirname:
           The directory where the octopus data was produced

        iteration:
           Iteration to make the difference, if absent the last iteration
           will be taken
        spin:
           Integer with the spin for which the data will be taken, if absent
           the addition of both spins will be returned

    Returns:
        A numpy array for the number of atoms with the Hirsfeld charges

    """
    if iteration is None:
        iteration = get_last_iteration(dirname)

    inp = pychemia.code.octopus.OctopusInput(dirname + '/inp')
    if inp.variables['SpinComponents'] == 'spin_polarized':
        if spin is None:
            it, hirshfeld1 = hirshfeld(dirname, iteration, spin=1)
            it, hirshfeld2 = hirshfeld(dirname, iteration, spin=2)
            return it, hirshfeld1 + hirshfeld2
        else:
            sp = '-sp' + str(spin)
    else:
        sp = ''

    filename = dirname + '/td.' + str(iteration).zfill(7) + '/density' + sp + '-hirschfeld.dat'
    if os.path.isfile(filename):
        hirshfeld_value = loadtxt(filename)
    else:
        filename = dirname + '/td.' + str(iteration).zfill(7) + '/density' + sp + '-hirshfeld.dat'
        hirshfeld_value = loadtxt(filename)

    return iteration, hirshfeld_value


def deflection(dirname, iteration=None):
    """
    Get the deflection angle respect to the initial velocity
    for all the atoms and for a given iteration and spin.

    Args:
        dirname:
           The directory where the octopus data was produced

        iteration:
           Iteration to make the difference, if absent the last iteration
           will be taken

    Returns:
        A numpy array for the number of atoms with the Hirsfeld charges

    """

    if iteration is None:
        iteration = get_last_iteration(dirname)

    natom, iterations, times, positions, velocities, forces = pychemia.code.octopus.coordinates(dirname)
    angle = zeros(natom)
    univel = zeros((natom, 3))

    if iteration == 0:
        print('Error: To get the deflection some evolution is necessary (iteration!=0)')

    else:

        # Vector of displacement from (0 -> t)
        dis = positions[iteration] - positions[0]
        # Magnitude of the displacement
        mag = apply_along_axis(linalg.norm, 1, dis)
        # Unitary vectors of displacement
        unidis = dis / mag[:, newaxis]

        # Unitary vector of velocity
        vel = array(velocities[0])
        for i in range(natom):
            # Magnitude of the velocities
            mag[i] = linalg.norm(vel[i])
            # Unitary vectors of velocity
            if mag[i] != 0:
                univel[i] = (1.0 / mag[i]) * vel[i]
                angle[i] = acos(dot(univel[i], unidis[i]))
            else:
                angle[i] = 0.0

    return angle


def magnitude_velocity(dirname, iteration=None):
    """
    Magnitude of velocity for each atom

    :param dirname: (str) Directory to analyse velocities
    :param iteration: (int) Number of iteration to use
    :return:
    """

    if iteration is None:
        iteration = get_last_iteration(dirname)

    natom, iterations, times, positions, velocities, forces = pychemia.code.octopus.coordinates(dirname)

    magvel = apply_along_axis(linalg.norm, 1, velocities[iteration])

    return magvel


def change_kinetic_energy(dirname, iteration=None):
    """
    Change of Kinetic energy respect to initial value

    :param dirname: (str) Directory to read
    :param iteration: (int) Iteration number to read
    :return:
    """
    if iteration is None:
        iteration = get_last_iteration(dirname)

    iterations, times, energy_dict = pychemia.code.octopus.energy(dirname)

    return energy_dict['Kinetic'][iteration] - energy_dict['Kinetic'][0]


def value_in_sphere(dirname, keys, iteration=None, radius=4.5, spin=None):
    """
    Compute the value of a certain scalar quantity inside a sphere of a
    given radius

    :param dirname: (str) Directory to read
    :param keys: (list) Variables to read the scalar quantity
    :param iteration: (int) Number of iteration to read
    :param radius: (float) Radius of the sphere to integrate
    :param spin: (int) Spin contribution to read
    :return:
    """
    if iteration is None:
        iteration = get_last_iteration(dirname)

    natom, iterations, times, positions, velocities, forces = pychemia.code.octopus.coordinates(dirname)

    # Is spin polarized?
    inp = pychemia.code.octopus.OctopusInput(dirname + '/inp')
    if inp.variables['SpinComponents'] == 'spin_polarized':
        if spin is None:
            it, cis1 = value_in_sphere(dirname, keys, iteration, radius, spin=1)
            it, cis2 = value_in_sphere(dirname, keys, iteration, radius, spin=2)
            return it, cis1 + cis2
        else:
            sp = '-sp' + str(spin)
    else:
        sp = ''

    cis = zeros(natom)
    for i in range(natom):

        filename = dirname + '/td.' + str(iteration).zfill(7)
        filename += '/density' + sp + '-analysis' + str(i + 1).zfill(3) + '.nc'

        if os.path.isfile(filename):
            data = scipy.io.netcdf.netcdf_file(filename, mmap=False)
            if (not keys[0] in data.variables) or (not keys[1] in data.variables):
                print('ERROR with ', keys)
            x = data.variables[keys[0]][:]
            if len(x) == 0:
                print('ERROR File empty:', filename)
            index = (abs(x - radius)).argmin()
            # print 'dirname',dirname
            # print 'index=',index
            # print 'radius (min,max)',X[0],X[-1]
            # Polynomial fit
            extra = 100
            if index + extra > len(x):
                final = len(x) - 1
            else:
                final = index + extra
            x = data.variables[keys[0]][index - extra:final]
            y = data.variables[keys[1]][index - extra:final]
            try:
                z = polyfit(x, y, 3)
                p = poly1d(z)
            except ValueError:
                print('Error fitting value for:', filename)

            # Get fit
            cis[i] = p(radius)
            data.close()
        else:
            print('ERROR Missing file:', filename)

    return iteration, cis


def charge_in_sphere(dirname, iteration=None, radius=4.5, spin=None):
    """
    Get the value at a given index of the
    data stored in the NetCDF file

    :param dirname: (str) Directory to read
    :param iteration: (int) Number of iteration to read
    :param radius: (float) Radius of the sphere to integrate
    :param spin: (int) Spin contribution to read
    :return:

    """

    keys = ['radeff', 'densit']

    return value_in_sphere(dirname, keys, iteration, radius, spin)


def dipolx_in_sphere(dirname, iteration=None, radius=4.5, spin=None):
    """
    Get the value at a given index of the
    data stored in the NetCDF file

    :param dirname: (str) Directory to read
    :param iteration: (int) Number of iteration to read
    :param radius: (float) Radius of the sphere to integrate
    :param spin: (int) Spin contribution to read
    :return:

    """

    keys = ['radeff', 'dipolX']

    return value_in_sphere(dirname, keys, iteration, radius, spin)
