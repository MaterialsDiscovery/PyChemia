from .searcher import Searcher
from pychemia import pcm_log


class GreyWolf(Searcher):
    def __init__(self, population, params=None, generation_size=32, stabilization_limit=10):

        Searcher.__init__(self, population, generation_size, stabilization_limit)
        self.a = None
        self.c = None
        self.set_params(params)

    def set_params(self, params):
        self.a = 1
        self.c = 1
        if params is None:
            params = {}
        if 'a' in params:
            self.a = params['a']
        if 'c' in params:
            self.c = params['c']

    def get_params(self):
        return {'a': self.a, 'c': self.c}

    def run_one(self):
        # Get a static selection of the values in the generation that are relaxed
        selection = self.population.ids_sorted(self.population.actives_evaluated)
        pcm_log.info('Size of selection : %d' % len(selection))
