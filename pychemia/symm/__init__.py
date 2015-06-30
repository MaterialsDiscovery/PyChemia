"""
Routines related to symmetry identification and manipulation.
The routines and classes heavily rely on spglib
"""

# __all__ = filter(lambda s: not s.startswith('_'), dir())

from pychemia import Structure

try:
    import pyspglib
    try:
        import pyspglib.spglib
        from pychemia.symm.symmetry import StructureSymmetry

        # Testing version of spglib
        st = Structure(symbols=['H'])
        symm = StructureSymmetry(st)
        ret = pyspglib.spglib.spg.dataset(symm._transposed_cell, symm._reduced, symm._numbers, 1e-5, -1.0)
        if type(ret[3]) is list:
            USE_SPGLIB = False
            print 'SPGLIB is present but outdated, please install spglib > 1.7'
        else:
            USE_SPGLIB = True
    except ImportError:
        print 'SPGLIB not found, symmetry module disabled'
        USE_SPGLIB = False

except ImportError:
    print 'SPGLIB not found, symmetry module disabled'
    USE_SPGLIB = False


if USE_SPGLIB:
    from .symmetry import StructureSymmetry, symmetrize
