from .searcher import Searcher


class MinimaHoppingMethod(Searcher):
    """
        class MinimaHopping
    """

    def __init__(self, population, params, generation_size=32, stabilization_limit=10):
        self.population = population
        Searcher.__init__(self, self.population, generation_size, stabilization_limit)
        self.params = None
        self.set_params(params)

    def set_params(self, params):
        self.params = params

    def run_one(self):
        pass

    def get_params(self):
        pass
