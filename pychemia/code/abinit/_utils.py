#!/usr/bin/env python

"""
This module provides general routines used by abipython
but not requiring the abipython classes
"""
from _htmlparser import MyHTMLParser

__author__ = "Guillermo Avendano-Franco"
__copyright__ = "Copyright 2012"
__version__ = "1.1"
__maintainer__ = "Guillermo Avendano-Franco"
__email__ = "guillermo.avendano@uclouvain.be"
__status__ = "Development"
__date__ = "Aug 27, 2012"

import os as _os
from ftplib import FTP as _FTP
from scipy.io import netcdf_file as _netcdf_file
from numpy import array as _array

from pychemia.utils.constants import angstrom_bohr as _angstrom_bohr
from pychemia.utils.periodic import atomic_symbol as _atom_symbol
from _input import InputVariables


def netcdf2dict(filename):
    """
    Read a NetCDF file and create a python dictionary with
    numpy arrays for variables

    Args:
        filename:
            NetCDF filename
    """
    if not _os.path.isfile(filename):
        print('ERROR: No such file: ', filename)
        return None
    output = {}
    netcdffile = _netcdf_file(filename, 'r', mmap=False)
    for ii in netcdffile.variables.keys():
        output[ii] = netcdffile.variables[ii][:]
    netcdffile.close()
    return output


def psp_name(atomicnumber, exchange, kind):
    """
    Return the filename of a certain PSP given the
    atomic number, exchange ('LDA' or 'GGA')
    and kind ('FHI','TM')
    """
    atom_symbol = str(_atom_symbol(atomicnumber))
    if kind == 'FHI' and exchange == 'LDA':
        filename = str(atomicnumber).zfill(2) + '-' + atom_symbol + '.LDA.fhi'
    elif kind == 'FHI' and exchange == 'GGA':
        filename = str(atomicnumber).zfill(2) + '-' + atom_symbol + '.GGA.fhi'
    elif kind == 'TM' and exchange == 'LDA':
        filename = str(atomicnumber) + atom_symbol.lower() + '.pspnc'
    elif kind == 'DEN' and exchange == 'AE':
        filename = '0.' + str(atomicnumber).zfill(2) + '-' + atom_symbol + '.8.density.AE'
    else:
        print('Not know kind of PSP')
        filename = ''
    return filename


def get_ftp_psp(atomicnumber, exchange, kind):
    """
    Get the ftp path to get the PSP
    """
    if kind == 'FHI' and exchange == 'LDA':
        ftppath = '/pub/abinitio/Psps/LDA_FHI/'
    elif kind == 'FHI' and exchange == 'GGA':
        ftppath = '/pub/abinitio/Psps/GGA_FHI/'
    elif kind == 'TM' and exchange == 'LDA':
        ftppath = '/pub/abinitio/Psps/LDA_TM.psps/' + str(atomicnumber).zfill(2) + '/'
    elif kind == 'DEN' and exchange == 'AE':
        ftppath = '/pub/abinitio/Psps/AE_DEN/'
    else:
        print('Not know kind of PSP')
        ftppath = ''
    return ftppath


def get_all_psps(basedir, exchange, kind):
    directory = basedir + '/' + exchange + '_' + kind
    if not _os.path.isdir(directory):
        _os.mkdir(directory)

    ftp = _FTP('ftp.abinit.org')  # connect to host, default port
    ftp.login()  # user anonymous, passwd anonymous@
    for i in range(1, 113):
        filename = psp_name(i, exchange, kind)
        if not _os.path.isfile(directory + '/' + filename):
            print('Getting...' + filename)
            nofile = True
            while nofile:
                try:
                    res = ftp.retrbinary('RETR ' + get_ftp_psp(i, exchange, kind) + filename,
                                         open(directory + '/' + filename, 'wb').write)
                    if _os.path.getsize(directory + '/' + filename) == 0:
                        _os.remove(directory + '/' + filename)
                        nofile = False
                    else:
                        nofile = False
                except ValueError:
                    print('Could not download ' + filename)
                    ftp.close()
                    if _os.path.isfile(directory + '/' + filename):
                        _os.remove(directory + '/' + filename)
                    ftp = _FTP('ftp.abinit.org')  # connect to host, default port
                    ftp.login()  # user anonymous, passwd anonymous@
                    nofile = False
    ftp.close()


def xyz2input(filename):
    """
    Reads a .xyz and return an ABINIT input
    as a python dictionary
    """

    abiinput = InputVariables()
    atomdict = _atom_symbol()
    rf = open(filename, 'r')

    natom = int(rf.readline())
    typat = []
    znucl = []
    xcart = []
    xangst = []

    ntypat = 0
    rf.readline()
    data = rf.readlines()
    for i in range(natom):
        atom = data[i].split()
        atomnumber = atomdict[atom[0]]
        if atomnumber not in znucl:
            ntypat += 1
            znucl.append(atomnumber)
        typat.append(znucl.index(atomnumber) + 1)
        xangst = xangst + [float(atom[1]), float(atom[2]), float(atom[3])]

    abiinput.variables['natom'] = _array([natom])
    abiinput.variables['znucl'] = _array(znucl)
    abiinput.variables['ntypat'] = _array([ntypat])
    abiinput.variables['typat'] = _array(typat)
    abiinput.variables['xcart'] = _angstrom_bohr * _array(xangst)

    return abiinput


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


def merge(abi_into, abi_from, filename=None):
    """
    From the variables present in the file 'abi_into'
    will try to recover the value from 'abi_from'
    if the value has changed, the value will be updated.

    The new variables will be save in the given filename
    If filename==None the new values will overwrite
    abi_into
    
    Example:
    
       merge('abinit.in','abinit_xo_OUT.nc')

    It will update the abinit input with the values of the output

    """
    abinit_into = InputVariables(abi_into)
    abinit_from = InputVariables(abi_from)

    for i in abinit_into.variables.keys():
        if i in abinit_from.variables.keys():
            abinit_into.variables[i] = abinit_from.variables[i]

    if filename is None:
        filename = abi_into

    abinit_into.write(filename)


def plot_simple(variables, varname):
    from matplotlib.pylab import subplots
    from numpy import arange, mean, apply_along_axis, linalg
    from math import sqrt

    fig, ax = subplots(nrows=1, ncols=1)
    fig.set_size_inches(15, 4)

    ndtset = variables['ndtset'][0]
    lens = _array([len(variables[varname + str(x + 1)]) for x in range(ndtset)])

    x = arange(ndtset) + 1

    if max(lens) == min(lens):
        if lens[0] == 1:
            y = _array([variables['etotal' + str(x + 1)][0] for x in range(ndtset)])
        elif lens[0] % 3 == 0:
            y = _array([mean(apply_along_axis(linalg.norm, 1, variables['fcart' + str(x + 1)].reshape((-1, 3))))
                        for x in range(ndtset)])
        else:
            y = _array([sqrt(sum(variables['fcart' + str(x + 1)] ** 2)) for x in range(ndtset)])
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
