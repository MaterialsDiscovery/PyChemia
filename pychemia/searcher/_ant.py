__author__ = 'guilleaf'

from _genealogy import Genealogy


class AntColony(Genealogy):

    def __init__(self, population, evaluator, objective_function, stabilization_limit=10, fraction_evaluated=0.8,
                 param=1):
        self.population = population
        self.evaluator = evaluator
        self.objective_function = objective_function
        self.stabilization_limit = stabilization_limit
        self.fraction_evaluated = fraction_evaluated
        Genealogy.__init__(self, self.population, self.evaluator)
        self.param = None
        self.set_param(param)

    def set_param(self, param):
        self.param = param

    def run_one_cycle(self):
        pass

    def run_all_cycles(self):
        pass
