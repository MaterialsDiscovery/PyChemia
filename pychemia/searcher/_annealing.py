
from _genealogy import Genealogy
from _searcher import Searcher

class SimulatedAnnealing(Genealogy, Searcher):

    def __init__(self, population, evaluator, objective_function, params, stabilization_limit=10,
                 fraction_evaluated=0.8):
        self.population = population
        self.evaluator = evaluator
        self.objective_function = objective_function
        self.stabilization_limit = stabilization_limit
        self.fraction_evaluated = fraction_evaluated
        Genealogy.__init__(self, self.population, self.evaluator)
        self.params = None
        self.set_params(params)

    def set_params(self, params):
        self.params = params

    def run_one_cycle(self):
        pass

    def run_all_cycles(self):
        pass
