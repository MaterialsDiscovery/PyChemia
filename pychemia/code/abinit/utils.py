import os
import numpy as np
from scipy.io import netcdf_file
from pychemia.utils.periodic import atomic_symbol
from .htmlparser import MyHTMLParser

"""
This module provides general routines used by abipython
but not requiring the abipython classes
"""

__author__ = "Guillermo Avendano-Franco"
__copyright__ = "Copyright 2016"
__version__ = "1.1"
__maintainer__ = "Guillermo Avendano-Franco"
__email__ = "guillermo.avendano@uclouvain.be"
__status__ = "Development"
__date__ = "May 13, 2016"


def netcdf2dict(filename):
    """
    Read a NetCDF file and create a python dictionary with
    numpy arrays for variables

    Args:
        filename:
            NetCDF filename
    """
    if not os.path.isfile(filename):
        print('ERROR: No such file: ', filename)
        return None
    output = {}
    netcdffile = netcdf_file(filename, 'r', mmap=False)
    for ii in netcdffile.variables.keys():
        output[ii] = netcdffile.variables[ii][:]
    netcdffile.close()
    return output


def psp_name(atomicnumber, exchange, kind):
    """
    Return the filename of a certain PSP given the
    atomic number, exchange ('LDA' or 'GGA')
    and kind ('FHI','TM')

    :param atomicnumber: (int) Atomic number
    :param exchange: (str) 'LDA' or 'GGA'
    :param kind: (str) Source of Pseudopotentials
    :return:
    """
    atom_symbol = str(atomic_symbol(atomicnumber))
    if kind == 'FHI' and exchange == 'LDA':
        filename = str(atomicnumber).zfill(2) + '-' + atom_symbol + '.LDA.fhi'
    elif kind == 'FHI' and exchange == 'GGA':
        filename = str(atomicnumber).zfill(2) + '-' + atom_symbol + '.GGA.fhi'
    elif kind == 'CORE' and exchange == 'LDA':
        filename = str(atomicnumber) + atom_symbol.lower() + '.1s_psp.mod'
    elif kind == 'GTH' and exchange == 'LDA':
        filename = str(atomicnumber).zfill(2) + atom_symbol.lower() + '.pspgth'
    elif kind == 'TM' and exchange == 'LDA':
        filename = str(atomicnumber) + atom_symbol.lower() + '.pspnc'
    elif kind == 'AE' and exchange == 'DEN':
        filename = '0.' + str(atomicnumber).zfill(2) + '-' + atom_symbol + '.8.density.AE'
    elif kind == 'FC' and exchange == 'DEN':
        filename = str(atomicnumber).zfill(2) + '-' + atom_symbol + '.8.fc'
    elif kind == 'PAW' and exchange == 'GGA':
        filename = 'JTH-PBE-atomicdata-0.2/' + atom_symbol + '.GGA_PBE-JTH.xml'
    elif kind == 'PAW' and exchange == 'LDA':
        filename = 'JTH-LDA-atomicdata-0.2/' + atom_symbol + '.LDA_PW-JTH.xml'
    elif kind == 'HGH' and exchange == 'GGA':
        filename = str(atomicnumber).zfill(2) + atom_symbol.lower() + '.pbe_hgh'
    elif kind == 'ONC' and exchange == 'PBE':
        filename = 'pbe_s_sr' + os.sep + atom_symbol + '.psp8'
    else:
        print('The combination of exchange=%s and kind=%s is not known' % (exchange, kind))
        filename = ''
    return filename


def split_varname(varname):
    if varname[-2:].isdigit():
        prefix = varname[:-2]
        suffix = varname[-2:]
    elif varname[-2].isdigit() and varname[-1] == '?':
        prefix = varname[:-2]
        suffix = varname[-2:]
    elif varname[-1].isdigit() and varname[-2] == '?':
        prefix = varname[:-2]
        suffix = varname[-2:]
    elif varname[-1].isdigit():
        prefix = varname[:-1]
        suffix = varname[-1:]
    elif varname[-1] == '?':
        prefix = varname[:-1]
        suffix = varname[-1:]
    else:
        prefix = varname
        suffix = ''
    return prefix, suffix


def plot_simple(variables, varname):
    from matplotlib.pylab import subplots
    from numpy import arange, mean, apply_along_axis, linalg
    from math import sqrt

    fig, ax = subplots(nrows=1, ncols=1)
    fig.set_size_inches(15, 4)

    ndtset = variables['ndtset'][0]
    lens = np.array([len(variables[varname + str(x + 1)]) for x in range(ndtset)])

    x = arange(ndtset) + 1

    if max(lens) == min(lens):
        if lens[0] == 1:
            y = np.array([variables['etotal' + str(x + 1)][0] for x in range(ndtset)])
        elif lens[0] % 3 == 0:
            y = np.array([mean(apply_along_axis(linalg.norm, 1, variables['fcart' + str(x + 1)].reshape((-1, 3))))
                          for x in range(ndtset)])
        else:
            y = np.array([sqrt(sum(variables['fcart' + str(x + 1)] ** 2)) for x in range(ndtset)])

        ax.plot(x, y, 'ro')
        ax.plot(x, y, 'b-')
        ax.set_xlabel('Dataset')
        ax.set_ylabel(varname)
        ax.set_xlim(1, ndtset + 1)
        if ndtset < 30:
            ax.set_xticks(arange(ndtset + 1))
            ax.grid(which='major', axis='both')


def abihelp(varname):
    import json

    hp = MyHTMLParser()

    import pychemia.code.abinit as _pdca

    rf = open(_pdca.__path__[0] + '/ABINIT_variables.json', 'r')
    variables = json.load(rf)
    rf.close()

    if varname not in variables.keys():
        print('ERROR: ' + varname + '  is not in the list of variables of ABINIT')
        return
    else:
        abivar = variables[varname]
        print(varname)
        print('')
        print('DEFINITION:', abivar['definition'])
        print('SECTION:   ', abivar['section'])
        print('DEFAULT:   ', hp.feed(abivar['default']))
        print('')
        print(hp.feed(abivar['text']))
        print('')
