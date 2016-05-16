"""
PyChemia works with atomic structures, this module provides two main classes to store structure information.
The class 'Structure' and the class 'Composition' are widely used for the entire package.

"""

from .composition import Composition
from .structure import Structure


# __all__ = filter(lambda s: not s.startswith('_'), dir())
