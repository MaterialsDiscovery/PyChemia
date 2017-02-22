"""
The classes and methods of this package, take one or several structures and based only on geometrical
properties, computes fingerprints, modify structures, match atoms from two structures, split and create
surfaces.
The functionality provided by this package is critical for many of the operations done by global search
methods applied for structural search.
"""

from .analysis import StructureAnalysis
from .changer import StructureChanger
from .matching import StructureMatch
from .cluster import ClusterAnalysis, ClusterMatch
from .surface import rotate_along_indices
from . import splitting
