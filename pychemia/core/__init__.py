"""
pychemia.core is the central subpackage in pychemia. It offers three important classes that are widely
used by all other modules.

The class 'Composition' stores and manipulates chemical compositions.

The class 'Structure' is used to store and manipulated chemical structures, ie, a composition with associated geometry
and in the case of crystals, a lattice.

The class 'Element' allow the retrieval of chemical and physical information about chemical species.
"""

from .composition import Composition
from .structure import Structure
from .element import Element

# __all__ = filter(lambda s: not s.startswith('_'), dir())
