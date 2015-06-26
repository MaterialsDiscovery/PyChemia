"""
PyChemia is an open-source Python Library for materials structural search. The purpose of the code is to create
a method agnostic framework for materials discovery and design using a variety of methods from Metaheuristic to
Dynamical such as minima hoping method (MHM)
"""

import logging
pcm_log = logging.getLogger(__name__)
pcm_log.addHandler(logging.NullHandler())

from core import Structure, Composition, Lattice
import analysis
import calc
import code
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
import population

from _info import __author__, __copyright__, __version__, __email__, __status__, __date__, version

# __all__ = filter(lambda s: not s.startswith('_'), dir())

def info():

    print 'Pychemia\n--------\n'
    print 'Version: '+__version__
    print 'Path:    '+__package__
    print 'Date:    '+__date__
    print

    import sys
    print 'Python version='+sys.version+'\n'

    try:
        mm = __import__('pymongo')
        print '%10s %10s   %s' % ('pymongo', mm.version, mm.__path__[0])
    except ImportError:
        print 'pymongo Not Found'

    for modui in ['numpy', 'scipy', 'mayavi', 'Scientific', 'matplotlib', 'pyhull', 'pymatgen', 'qmpy']:
        try:
            mm = __import__(modui)
            print '%10s %10s   %s' % (modui, mm.__version__, mm.__path__[0])
        except ImportError:
            print '%10s Not Found' % modui

    try:
        import ase
        from ase import version as ase_version
        print '%10s %10s   %s' % ('ase', ase_version.version_base, ase.__path__[0])
    except ImportError:
        print 'ase Not Found'
