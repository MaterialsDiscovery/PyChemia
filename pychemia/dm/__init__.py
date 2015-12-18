"""
Routines related to Data mining
"""
try:
    from network import NetworkAnalysis
except ImportError:
    print 'Networkx not present'

# __all__ = filter(lambda s: not s.startswith('_'), dir())
