"""
PyChemia is a python library for automatize atomistic-level calculations. The library provide an API to manipulate
structures, store structures and properties in a Mongo database, search for new structures of a given composition,
interact with several atomistic simulation codes and visualize atomistic-related data

"""

import os
import sys
import logging
import json
import importlib

with open(__path__[0] + os.sep + 'setup.json') as rf:
    data = json.load(rf)

__author__ = data['author']
__copyright__ = data['copyright']
__version__ = data['version']
__email__ = data['email']
__status__ = data['status']
__date__ = data['date']


spec = importlib.util.find_spec("scipy")
HAS_SCIPY = spec is not None



try:
    try:
        import spglib as spg
    except ImportError:
        from pyspglib import spglib as spg
    HAS_SPGLIB = True
except ImportError:
    HAS_SPGLIB = False

spec = importlib.util.find_spec("matplotlib")
HAS_MATPLOTLIB = spec is not None

if HAS_MATPLOTLIB and 'matplotlib' not in sys.modules:
    import matplotlib
    matplotlib.use('agg')

spec = importlib.util.find_spec("matplotlib")
HAS_MATPLOTLIB = spec is not None

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
