from _function import FunctionEvaluator, FunctionObjectiveFunction
from abc import ABCMeta, abstractmethod, abstractproperty


class Evaluator:
    __metaclass__ = ABCMeta

    def __init__(self):
        self.population = None

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
