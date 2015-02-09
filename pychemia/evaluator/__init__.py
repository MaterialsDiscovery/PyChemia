__author__ = 'Guillermo Avendano Franco'

#from _ase import AseEvaluator, AseObjectiveFunction
from _function import FunctionEvaluator, FunctionObjectiveFunction

from abc import ABCMeta, abstractmethod, abstractproperty


class Evaluator():
    __metaclass__ = ABCMeta

    def __init__(self):
        pass

    @abstractmethod
    def initialize(self, population):
        self.population = population

    @abstractmethod
    def run(self):
        pass

    @abstractmethod
    def stop(self):
        pass

    @abstractproperty
    def is_running(self):
        pass

class EvaluatorDummy(Evaluator):

    def __init__(self):
        pass

    def initialize(self, population):
        self.population = population

    def run(self):
        pass

    def stop(self):
        pass

    @property
    def is_running(self):
        return True
