"""
PyChemia is a python library for automatize atomistic-level calculations. The library provide an API to manipulate
structures, store structures and properties in a Mongo database, search for new structures of a given composition,
interact with several atomistic simulation codes and visualize atomistic-related data

"""

import os
import sys
import logging
import json

with open(__path__[0] + os.sep + 'setup.json') as rf:
    data = json.load(rf)

__author__ = data['author']
__copyright__ = data['copyright']
__version__ = data['version']
__email__ = data['email']
__status__ = data['status']
__date__ = data['date']

try:
    import scipy

    HAS_SCIPY = True
except ImportError:
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
    if 'matplotlib' not in sys.modules:
        matplotlib.use('agg')

except ImportError:
    HAS_MATPLOTLIB = False

try:
    import mayavi
    HAS_MAYAVI = True
except ImportError:
    HAS_MAYAVI = False

try:
    import vtk
    HAS_VTK = True
except ImportError:
    HAS_VTK = False

try:
    import pyhull

    HAS_PYHULL = True
except ImportError:
    HAS_PYHULL = False

try:
    import networkx

    HAS_NETWORKX = True
except ImportError:
    HAS_NETWORKX = False

try:
    import pymongo

    if pymongo.version_tuple[0] < 3:
        HAS_PYMONGO = False
    else:
        HAS_PYMONGO = True
except ImportError:
    pymongo = None
    HAS_PYMONGO = False

try:
    import gridfs

    HAS_GRIDFS = True
except ImportError:
    gridfs = None
    HAS_GRIDFS = False

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
    Show basic information about PyChemia, its location and version.
        Also information about other libraries used by PyChemia
        both mandatory and optional

    """
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

    for modui in ['numpy', 'scipy', 'spglib', 'matplotlib', 'nose', 'coverage', 'pyhull', 'pymatgen',
                  'networkx', 'ase', 'mayavi', 'qmpy', ]:
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


class Version:
    @staticmethod
    def full_version():
        return 'PyChemia version: ' + __version__ + ' from: ' + __date__

    def __init__(self):
        pass
