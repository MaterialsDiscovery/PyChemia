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

import os
import numpy as np
from ftplib import FTP as _FTP
from scipy.io import netcdf_file as _netcdf_file
from pychemia.utils.periodic import atomic_symbol


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
    atom_symbol = str(atomic_symbol(atomicnumber))
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
    if not os.path.isdir(directory):
        os.mkdir(directory)

    ftp = _FTP('ftp.abinit.org')  # connect to host, default port
    ftp.login()  # user anonymous, passwd anonymous@
    for i in range(1, 113):
        filename = psp_name(i, exchange, kind)
        if not os.path.isfile(directory + '/' + filename):
            print('Getting...' + filename)
            nofile = True
            while nofile:
                try:
                    res = ftp.retrbinary('RETR ' + get_ftp_psp(i, exchange, kind) + filename,
                                         open(directory + '/' + filename, 'wb').write)
                    if os.path.getsize(directory + '/' + filename) == 0:
                        os.remove(directory + '/' + filename)
                        nofile = False
                    else:
                        nofile = False
                except ValueError:
                    print('Could not download ' + filename)
                    ftp.close()
                    if os.path.isfile(directory + '/' + filename):
                        os.remove(directory + '/' + filename)
                    ftp = _FTP('ftp.abinit.org')  # connect to host, default port
                    ftp.login()  # user anonymous, passwd anonymous@
                    nofile = False
    ftp.close()


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
