"""
Routines related to Report generation
"""
from pychemia import HAS_MATPLOTLIB
from .structure_plot import StructurePlot

from .dos import DensityOfStates, plot_one_dos, plot_many_dos
# from _pyprocar import BandStructure

if HAS_MATPLOTLIB:
    from pychemia.visual.searcher import plot_evolution_circular

# __all__ = filter(lambda s: not s.startswith('_'), dir())
