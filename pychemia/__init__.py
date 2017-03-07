"""
PyChemia is a python library for automatize atomistic-level calculations. The library provide an API to manipulate
structures, store structures and properties in a Mongo database, search for new structures of a given composition,
interact with several atomistic simulation codes and visualize atomistic-related data

"""

from __future__ import print_function
import sys
import logging

__author__ = "Guillermo Avendano-Franco"
__copyright__ = "Copyright 2017"
__version__ = "0.17.3"
__email__ = "gufranco@mail.wvu.edu"
__status__ = "Development"
__date__ = "March 7, 2017"

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
    import Scientific
    HAS_SCIENTIFIC = True
except ImportError:
    HAS_SCIENTIFIC = False

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
from .core import Structure, Composition
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
from pychemia.crystal import samples


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

    for modui in ['numpy', 'scipy', 'mayavi', 'Scientific', 'matplotlib',
                  'future', 'nose', 'coverage', 'spglib', 'pyhull', 'pymatgen', 'qmpy', ]:
        try:
            mm = __import__(modui)
            print('%10s %10s   %s' % (modui, mm.__version__, mm.__path__[0]))
        except ImportError:
            print('%10s %10s Not Found' % (modui, ''))

    try:
        import ase
        from ase import version as ase_version
        print('%10s %10s   %s' % ('ase', ase_version.version_base, ase.__path__[0]))
    except ImportError:
        print('%10s %10s Not Found' % ('ase', ''))


class Version:
    @staticmethod
    def full_version():
        return 'PyChemia Version=' + __version__ + ' from=' + __date__

    def __init__(self):
        pass
