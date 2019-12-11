"""
PyChemia is a python library for automatize atomistic-level calculations. The library provides classes and routines to
manipulate the geometry of atomic structures, with emphasis on crystals.
It has connectors to manipulate the input for several DFT codes and interpret the output from those codes too.
PyChemia allows the store of structures and properties in a Mongo database, search for new structures of a given
composition, and perform a structural search using several metaheuristic global search algorithms included on PyChemia.
"""

import sys
import logging
import importlib

from pychemia.version import version as __version__
from pychemia.version import author as __author__
from pychemia.version import copyright as __copyright__
from pychemia.version import email as __email__
from pychemia.version import status as __status__
from pychemia.version import date as __date__

# Mandatory dependencies
########################

# These check have more sense during develompent
# Installing with pip will fulfill these dependencies

spec = importlib.util.find_spec("scipy")
HAS_SCIPY = spec is not None

# Versions 1.8.x or before used to be pyspglib
# Removed on (2019-12-11)
# try:
#     try:
#         import spglib as spg
#     except ImportError:
#         from pyspglib import spglib as spg
#     HAS_SPGLIB = True
# except ImportError:
#     HAS_SPGLIB = False

spec = importlib.util.find_spec("spglib")
HAS_SPGLIB = spec is not None

spec = importlib.util.find_spec("matplotlib")
HAS_MATPLOTLIB = spec is not None

if HAS_MATPLOTLIB and 'matplotlib' not in sys.modules:
    import matplotlib

    matplotlib.use('agg')

spec = importlib.util.find_spec("matplotlib")
HAS_MATPLOTLIB = spec is not None

spec = importlib.util.find_spec("psutil")
HAS_PSUTIL = spec is not None

# Optional Dependencies
#######################

spec = importlib.util.find_spec("mayavi")
HAS_MAYAVI = spec is not None

spec = importlib.util.find_spec("vtk")
HAS_VTK = spec is not None

spec = importlib.util.find_spec("pyhull")
HAS_PYHULL = spec is not None

spec = importlib.util.find_spec("networkx")
HAS_NETWORKX = spec is not None

spec = importlib.util.find_spec("pymongo")
HAS_PYMONGO = spec is not None

if HAS_PYMONGO:
    import pymongo

    if pymongo.version_tuple[0] < 3:
        HAS_PYMONGO = False
    else:
        HAS_PYMONGO = True

spec = importlib.util.find_spec("gridfs")
HAS_GRIDFS = spec is not None

spec = importlib.util.find_spec("ase")
HAS_ASE = spec is not None

spec = importlib.util.find_spec("pymatgen")
HAS_PYMATGEN = spec is not None

pcm_log = logging.getLogger(__name__)
pcm_log.addHandler(logging.NullHandler())

from .core import Structure, Composition, Element
from . import analysis
from . import db
from . import crystal
from . import io
from . import runner
from . import searcher
from . import utils
from . import web
from . import code
from . import population
from .core.from_file import structure_from_file


def info():
    """
    Show basic information about PyChemia, the version, location of the package and registered date of release
    Both mandatory and optional dependecies are shown, including version and location.

    """
    print('PyChemia\n--------\n')
    print('Version: ' + __version__)
    print('Path:    ' + __path__[0])
    print('Date:    ' + __date__)
    print()

    import sys
    print('Python version=' + sys.version + '\n')

    for modui in ['numpy', 'scipy', 'spglib', 'matplotlib', 'pymongo', 'psutil',
                  'nose', 'coverage', 'pyhull', 'pymatgen',
                  'networkx', 'ase', 'mayavi', 'qmpy', ]:
        if modui == 'numpy':
            print("Mandatory dependencies:")
        if modui == 'nose':
            print("\nOptional dependencies:")
        try:
            mm = __import__(modui)
            print('%10s %10s   %s' % (modui, mm.__version__, mm.__path__[0]))
        except ImportError:
            print('%10s %10s Not Found' % (modui, ''))

    try:
        import vtk
        print('%10s %10s   %s' % ('vtk', vtk.VTK_VERSION, vtk.__path__[0]))
    except ImportError:
        print('%10s %10s Not Found' % ('vtk', ''))
