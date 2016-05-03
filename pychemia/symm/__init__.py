"""
Routines related to symmetry identification and manipulation.
The routines and classes heavily rely on spglib
"""

# __all__ = filter(lambda s: not s.startswith('_'), dir())

from pychemia import Structure

try:

    try:
        import spglib as spg
    except ImportError:
        from pyspglib import spglib as spg

    from pychemia.symm.symmetry import StructureSymmetry

    # Testing version of spglib
    st = Structure(symbols=['H'])
    symm = StructureSymmetry(st)
    ret = spg.spglib.spg.dataset(symm.transposed, symm.reduced, symm.numbers, 1e-5, -1.0)
    if type(ret[3]) is list:
        USE_SPGLIB = False
        version = "%d.%d.%d" % spg.get_version()
        print('SPGLIB current version is %s, please install spglib > 1.9' % version)
    else:
        USE_SPGLIB = True

except ImportError:
    USE_SPGLIB = False

if USE_SPGLIB:
    from .symmetry import StructureSymmetry, symmetrize
