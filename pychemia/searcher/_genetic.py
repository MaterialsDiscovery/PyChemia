from _searcher import Searcher
from pychemia import pcm_log

"""
Implementation of Genetic Algorithms
This is the abstract layer
No references to Structural Search should be appear
This module should be general enough even for stock market prediction!
"""

__author__ = 'Guillermo Avendano-Franco'


class GeneticAlgorithm(Searcher):

    def __init__(self, population, params=None, fraction_evaluated=0.95, generation_size=32, stabilization_limit=10):
        Searcher.__init__(self, population, fraction_evaluated, generation_size, stabilization_limit)
        self.nelite = None
        self.mutation_prob = None
        self.cross_prob = None
        self.set_params(params)

    def set_params(self, params):
        self.nelite = 2
        self.mutation_prob = 0.5
        self.cross_prob = 0.5
        if params is None:
            params = {}
        if 'nelite' in params:
            self.nelite = params['nelite']
        if 'mutation_prob' in params:
            self.mutation_prob = params['mutation_prob']
        if 'cross_prob' in params:
            self.cross_prob = params['cross_prob']

    def get_params(self):
        return {'nelite': self.nelite, 'mutation_prob': self.mutation_prob, 'cross_prob': self.cross_prob}

    def run_one(self):
        # Get a static selection of the values in the generation that are relaxed
        selection = self.population.ids_sorted(self.population.actives_evaluated)
        pcm_log.info('Size of selection : %d' % len(selection))
