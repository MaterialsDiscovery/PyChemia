"""
Routines related to convert pymatgen structure into pymatdis structures
"""

from pychemia import HAS_PYMATGEN

if HAS_PYMATGEN:
    from ._pymatgen import cif2structure, pychemia2pymatgen, pymatgen2pychemia

# __all__ = filter(lambda s: not s.startswith('_'), dir())
