"""
Routines related to Report generation
"""
from pychemia import HAS_MATPLOTLIB, HAS_MAYAVI
from .dos import DensityOfStates, plot_one_dos, plot_many_dos
from .povray import StructurePovray

if HAS_MAYAVI:
    from .structure_plot import StructurePlot

# from _pyprocar import BandStructure

if HAS_MATPLOTLIB:
    from pychemia.visual.searcher import plot_evolution_circular

# __all__ = filter(lambda s: not s.startswith('_'), dir())
