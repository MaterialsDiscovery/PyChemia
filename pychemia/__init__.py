"""
PyChemia is an open-source Python Library for materials structural search. The purpose of the code is to create
a method agnostic framework for materials discovery and design using a variety of methods from Metaheuristic to
Dynamical such as minima hoping method (MHM)
"""
import sys

if 'matplotlib' not in sys.modules:
    import matplotlib

    matplotlib.use('agg')

import logging

pcm_log = logging.getLogger(__name__)
pcm_log.addHandler(logging.NullHandler())
from .core import Structure, Composition, Lattice
import analysis
import db
import dft
import dm
import gui
import io
import report
import runner
import searcher
import symm
import utils
import web
import external
import serializer
from _info import __author__, __copyright__, __version__, __email__, __status__, __date__, Version

try:
    import pymongo

    HAS_PYMONGO = True
except ImportError:
    HAS_PYMONGO = False

try:
    import scipy

    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False

try:
    import Scientific

    HAS_SCIENTIFIC = True
except ImportError:
    HAS_SCIENTIFIC = False

try:
    import pyspglib

    HAS_SPGLIB = True
except ImportError:
    HAS_SPGLIB = False

try:
    import ase

    HAS_ASE = True
except ImportError:
    HAS_ASE = False

try:
    import pymatgen

    HAS_PYMATGEN = True
except ImportError:
    HAS_PYMATGEN = False

try:
    import gridfs

    HAS_GRIDFS = True
except ImportError:
    HAS_GRIDFS = False

import code
import population
from .core.from_file import structure_from_file


def info():
    print('PyChemia\n--------\n')
    print('Version: ' + __version__)
    print('Path:    ' + __path__[0])
    print('Date:    ' + __date__)
    print()

    import sys
    print('Python version=' + sys.version + '\n')

    try:
        mm = __import__('pymongo')
        print('%10s %10s   %s' % ('pymongo', mm.version, mm.__path__[0]))
    except ImportError:
        print('pymongo Not Found')

    for modui in ['numpy', 'scipy', 'mayavi', 'Scientific', 'matplotlib', 'pyhull', 'pymatgen', 'qmpy']:
        try:
            mm = __import__(modui)
            print('%10s %10s   %s' % (modui, mm.__version__, mm.__path__[0]))
        except ImportError:
            print('%10s Not Found' % modui)

    try:
        import ase
        from ase import version as ase_version
        print('%10s %10s   %s' % ('ase', ase_version.version_base, ase.__path__[0]))
    except ImportError:
        print('ase Not Found')
