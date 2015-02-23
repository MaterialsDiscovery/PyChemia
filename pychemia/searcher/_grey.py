__author__ = 'Guillermo Avendano-Franco'

from _searcher import Searcher


class GreyWolf(Searcher):
    def __init__(self, population, params=None, fraction_evaluated=0.8, generation_size=32, stabilization_limit=10):
        self.population = population

        Searcher.__init__(self, self.population, fraction_evaluated, generation_size, stabilization_limit)
        self.a = None
        self.c = None
        self.set_params(params)

    def set_params(self, params):
        if params is None:
            self.a = 1
            self.c = 1
            return
        if 'a' not in params:
            self.a = 1
        else:
            self.a = params['a']
        if 'c' not in params:
            self.c = 1
        else:
            self.c = params['c']

    def run_one_cycle(self):
        pass
