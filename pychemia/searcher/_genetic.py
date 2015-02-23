"""
Implementation of Genetic Algorithms
This is the abstract layer
No references to Structural Search should be appear
This module should be general enough even for stock market prediction!
"""

__author__ = 'Guillermo Avendano-Franco'

from _searcher import Searcher


class GeneticAlgorithm(Searcher):

    def __init__(self, population, params, fraction_evaluated=0.8, generation_size=32, stabilization_limit=10):
        self.population = population
        Searcher.__init__(self, self.population, fraction_evaluated, generation_size, stabilization_limit)
        self.params = None
        self.set_params(params)

    def set_params(self, params):
        self.params = params

    def run_one_cycle(self):
        pass
