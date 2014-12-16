__author__ = 'viviane'

from abc import ABCMeta, abstractmethod


class Codes():
    __metaclass__ = ABCMeta

    def __init__(self):
        pass

    @abstractmethod
    def initialize(self, dirpath):
        pass

    @abstractmethod
    def set_inputs(self):
        pass

    @abstractmethod
    def get_ouputs(self):
        pass

    @abstractmethod
    def run(self):
        pass

    @abstractmethod
    def finalize(self):
        pass

    #@abstractproperty
    #def dirpath():
    #    pass
