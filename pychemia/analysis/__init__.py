"""
Set of classes that uses purely structural information to manipulate, split and match structures.
The classes and methods in this package are critical for many of the structural search functionality of PyChemia
"""

from .analysis import StructureAnalysis
from .changer import StructureChanger
from .matching import StructureMatch
from .cluster import ClusterAnalysis, ClusterMatch
from .surface import create_surface
from . import splitting
