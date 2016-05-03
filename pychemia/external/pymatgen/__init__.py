"""
Routines related to convert pymatgen structure into pymatdis structures
"""

try:
    from _pymatgen import cif2structure, pychemia2pymatgen, pymatgen2pychemia
except ImportError:
    print('The module pymatgen could not be imported')

# __all__ = filter(lambda s: not s.startswith('_'), dir())
