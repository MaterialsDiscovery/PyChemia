__author__ = 'Guillermo Avendano-Franco'

from pychemia.db import PyChemiaDB


class Population():
    """
    The class Population defines a set of structures stored in a Mongo Data Base
    This class is used by population-based optimization methods
    The population could be static or dynamic on their size.
    There are defined one or several leaders and the optimization function could be
    scalar or vector
    This class provides an abstraction to all those parameters and acts a an interface
    to manipulate the entries on the database.
    """

    def __init__(self, name):

        self.size = 0
        self._dynamic = True
        self._target_dim = 1
        self.db = PyChemiaDB(name)
        self.leaders = []
        self.active = []

    @property
    def is_static(self):
        return not self._dynamic

    @property
    def is_dynamic(self):
        return self._dynamic

    def add_leader(self, identifier):
        self.leaders.append(identifier)

    def del_leader(self, identifier):
        self.leaders.pop(identifier)

    def add_active(self, identifier):
        self.leaders.append(identifier)

    def del_active(self, identifier):
        self.leaders.pop(identifier)

    def replace_entry(self, identifier, structure):
        pass
