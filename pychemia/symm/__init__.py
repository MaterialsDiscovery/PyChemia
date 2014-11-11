"""
Routines related to Symmetry identification and manipulation
"""

#__all__ = filter(lambda s: not s.startswith('_'), dir())

try:
    import pyspglib._spglib as spg
    from _spglib import *
    USE_SPGLIB = True
except ImportError:
    print 'spglib functionality disabled'
    USE_SPGLIB = False
