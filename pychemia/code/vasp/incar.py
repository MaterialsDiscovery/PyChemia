import os
from .input import VaspInput


__author__ = "Guillermo Avendano-Franco"
__copyright__ = "Copyright 2016"
__version__ = "0.1"
__maintainer__ = "Guillermo Avendano-Franco"
__email__ = "gtux.gaf@gmail.com"
__status__ = "Development"
__date__ = "May 13, 2016"


def read_incar(filename='INCAR'):
    """
    Load the file INCAR in the directory 'path' or
    read directly the file 'path' and return an object
    'inputvars' for pychemia

    :param filename: (str) Filename of a INCAR file format
    :return:
    """
    if os.path.isfile(filename):
        filename = filename
    elif os.path.isdir(filename) and os.path.isfile(filename + '/INCAR'):
        filename += '/INCAR'
    else:
        raise ValueError('[ERROR] INCAR path not found: %s' % filename)

    iv = VaspInput(filename=filename)
    return iv


def write_incar(iv, filepath='INCAR'):
    """
    Takes an object inputvars from pychemia and
    save the file INCAR in the directory 'path' or
    save the file 'path' as a VASP INCAR file

    :param iv: (VaspInput) VASP Input variables
    :param filepath: (str) File path to write the INCAR file
    """

    if os.path.isdir(filepath):
        filename = filepath + '/INCAR'
    else:
        filename = filepath

    iv.write(filename)

