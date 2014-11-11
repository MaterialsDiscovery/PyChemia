"""
Population class
"""

__author__ = 'Guillermo Avendano-Franco'

import pychemia
if pychemia.db.USE_MONGO:
    from pychemia.db import PyChemiaDB
from pychemia import Structure, Composition


class Population():
    """
    The class Population defines a set of structures stored in a MongoDB database
    This class is used by population-based optimization methods
    The population could be static or dynamic on their size.
    There are defined one or several leaders and the optimization function could be
    scalar or vector
    This class provides an abstraction to all those parameters and acts a an interface
    to manipulate the entries on the database.
    """

    def __init__(self, name, new=False):
        """
        Init

        :param name:
        :param size:
        :param new:
        :return:
        """

        self.name = name
        self._dynamic = True
        self._target_dim = 1
        if pychemia.db.USE_MONGO:
            self.pcdb = PyChemiaDB(name)
        if new:
            self.pcdb.clean()
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

    def add_random(self, composition, size=1):
        """
        Add random structures to the population

        :param composition: (pychemia.Composition) The composition object
        :param size: (int) Number of element (default=1)
        """
        for i in range(size):
            structure = Structure().random_cell(composition)
            self.pcdb.insert(structure)

    def __repr__(self):
        return str(self.__class__)+'({0.name!r})'.format(self)

    def __str__(self):
        ret = 'Population\n\nName: {0.name!s}\n'.format(self)
        ret += 'Size: %d' % self.pcdb.entries.count()
        return ret

    def add_from_db(self, composition, dbname, sizemax=1):

        comp = Composition(composition)
        readdb = PyChemiaDB(dbname)

        index = 0
        for entry in readdb.entries.find({'formula': comp.formula, 'natom': comp.natom}):
            if index < sizemax:
                print 'Adding entry '+str(entry['_id'])+' from '+dbname
                self.pcdb.insert(Structure().fromdict(entry))
                index += 1

    def __len__(self):
        return self.pcdb.entries.count()

    @property
    def entries_ids(self):
        return [str(entry['_id']) for entry in self.pcdb.entries.find()]

    def set_all_active(self):

        for i in self.entries_ids:
            self.active.append(i)

    @property
    def actives(self):
        return self.active

    def update_entry(self, entry_id, structure):

        if pychemia.db.USE_MONGO:
            self.pcdb.update(entry_id, structure.todict())