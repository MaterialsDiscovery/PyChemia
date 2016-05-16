"""
Populations are objects collecting structures/properties usually for feeding global search algorithms.
Several populations are implemented.
*RelaxStructures* is used for structural search, *OrbitalDFTU* is used for optimization of orbital
 mixing and determination of U on DFT+U calculations. *NonCollinearMagMoms* optimize the directions
 of MagMoms for Non-Collinear-Spin calculations. *RealFunction* is a population of vector for optimize
 real-valuated functions. *LJCluster* stores LennardJones clusters.

"""

from pychemia import HAS_PYMONGO
from .realfunction import RealFunction

if HAS_PYMONGO:
    from .relaxstructures import RelaxStructures
    from .orbitaldftu import OrbitalDFTU, dmatpawu2params, params2dmatpawu, get_pattern, params_reshaped
    from .noncollinearmagmoms import NonCollinearMagMoms
    from .ljcluster import LJCluster
