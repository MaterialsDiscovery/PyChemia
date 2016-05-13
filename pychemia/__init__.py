"""
PyChemia is an open-source Python Library for materials structural search. The purpose of the code is to create
a method agnostic framework for materials discovery and design using a variety of methods from Metaheuristic to
Dynamical such as minima hoping method (MHM)
"""
from __future__ import print_function

try:
    import scipy

    HAS_SCIPY = True
except ImportError:
    # print("Library 'scipy' could not be found, several places of the code will be disabled")
    HAS_SCIPY = False

try:
    try:
        import spglib as spg
    except ImportError:
        from pyspglib import spglib as spg
    HAS_SPGLIB = True
except ImportError:
    HAS_SPGLIB = False


try:
    import matplotlib


    HAS_MATPLOTLIB = True
except ImportError:
    # print("Library 'matplotlib' could not be found, disabling visual functionality")
    HAS_MATPLOTLIB = False

try:
    import pyhull

    HAS_PYHULL = True
except ImportError:
    #print("Library 'pyhull' could not be found")
    HAS_PYHULL = False

try:
    import networkx

    HAS_NETWORKX = True
except ImportError:
    #print("Library 'networkx' could not be found, disabling pychemia.dm.NetworkAnalysis")
    HAS_NETWORKX = False

try:
    import Scientific
    HAS_SCIENTIFIC = True
except ImportError:
    #print("Library 'Scientific' could not be found")
    HAS_SCIENTIFIC = False

import sys

if HAS_MATPLOTLIB and 'matplotlib' not in sys.modules:
    import matplotlib
    matplotlib.use('agg')

try:
    import pymongo

    if pymongo.version_tuple[0] < 3:
        #print("Library 'pymongo' its too old, disabling pychemia.db.PyChemiaDB")
        HAS_PYMONGO = False
    else:
        HAS_PYMONGO = True
except ImportError:
    pymongo = None
    #print("Library 'pymongo' could not be found, disabling pychemia.db.PyChemiaDB")
    HAS_PYMONGO = False

try:
    import gridfs

    HAS_GRIDFS = True
except ImportError:
    gridfs = None
    #print("Library 'gridfs' could not be found, disabling pychemia.db.PyChemiaQueue")
    HAS_GRIDFS = False


try:
    import ase
    HAS_ASE = True
except ImportError:
    #print("Library 'ase' could not be found, disabling pychemia.external.ase")
    HAS_ASE = False

try:
    import pymatgen
    HAS_PYMATGEN = True
except ImportError:
    #print("Library 'pymatgen' could not be found, disabling pychemia.external.pymatgen")
    HAS_PYMATGEN = False

import logging

pcm_log = logging.getLogger(__name__)
pcm_log.addHandler(logging.NullHandler())
from .core import Structure, Composition
from . import analysis
from . import db
from . import crystal
from . import io
from . import runner
from . import searcher
from . import utils
from . import web
from ._info import __author__, __copyright__, __version__, __email__, __status__, __date__, Version
from . import code
from . import population
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
