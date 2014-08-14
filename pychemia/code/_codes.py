__author__ = 'viviane'

from abc import ABCMeta, abstractmethod

class Codes(ABC):
    __metaclass__ = ABCMeta

    @abstractmethod
    def initialize(self):
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

    @abstractproperty
    def dirpath(self):

