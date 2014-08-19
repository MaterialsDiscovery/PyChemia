__author__ = 'Guillermo Avendano-Franco'

from abc import ABCMeta


class MetaHeuristics():
    __metaclass__ = ABCMeta
    """
    Abstract class for all optimization algorithms that uses MetaHeuristics
    """

    def set_population(self, nsize):
        pass

