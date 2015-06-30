"""
Routines related to symmetry identification and manipulation.
The routines and classes heavily rely on spglib
"""

# __all__ = filter(lambda s: not s.startswith('_'), dir())

try:
    import pyspglib
    try:
        import pyspglib.spglib
        if 'primitive' not in pyspglib.spglib.__dict__ and 'primitive' in pyspglib.spglib.spg.__dict__:
            USE_SPGLIB = True
        else:
            print 'SPGLIB is present but outdated, please install spglib > 1.7'
            USE_SPGLIB = False
    except ImportError:
        print 'SPGLIB not found, symmetry module disabled'
        USE_SPGLIB = False

except ImportError:
    print 'SPGLIB not found, symmetry module disabled'
    USE_SPGLIB = False


if USE_SPGLIB:
    from symmetry import StructureSymmetry, symmetrize
