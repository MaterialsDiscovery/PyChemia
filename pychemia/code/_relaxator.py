__author__ = 'Guillermo Avendano-Franco'

from abc import ABCMeta, abstractmethod


class Relaxator():
    __metaclass__ = ABCMeta

    def __init__(self):
        pass

    @abstractmethod
    def set_params(self, params):
        pass

    @abstractmethod
    def run(self):
        pass

    @abstractmethod
    def get_forces_stress_energy(self):
        pass

    @abstractmethod
    def get_final_geometry(self):
        pass
