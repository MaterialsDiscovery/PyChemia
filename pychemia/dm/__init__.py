"""
Routines related to Data mining
"""
from pychemia import HAS_NETWORKX

if HAS_NETWORKX:
    from network import NetworkAnalysis

# __all__ = filter(lambda s: not s.startswith('_'), dir())
