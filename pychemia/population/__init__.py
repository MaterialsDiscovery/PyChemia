from pychemia import HAS_PYMONGO
from .realfunction import RealFunction

if HAS_PYMONGO:
    from .relaxstructures import RelaxStructures
    from .orbitaldftu import OrbitalDFTU
    from .noncollinearmagmoms import NonCollinearMagMoms
    from .ljcluster import LJCluster

"""
Populations are objects collecting structures/properties usually for feeding global search algorithms.
Several populations are implemented.

Several population classes are defined

__RelaxStructures__: Used for structural search.
__OrbitalDFTU__: Used for optimization of orbital mixing and determination of U on DFT+U calculations.
__NonCollinearMagMoms__: Optimize the directions of MagMoms for Non-Collinear-Spin calculations.
__RealFunction__: A population of vector for optimize real-valuated functions.
__LJCluster__: Stores Lennard-Jones clusters.

"""
