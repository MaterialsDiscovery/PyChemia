__author__ = 'Guillermo Avendano Franco'

from abc import ABCMeta, abstractmethod


class Searcher():

    __metaclass__ = ABCMeta

    @abstractmethod
    def set_params(self, params):
        pass

    @abstractmethod
    def run_one_cycle(self):
        pass