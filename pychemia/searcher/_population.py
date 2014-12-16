"""
Population class
"""

__author__ = 'Guillermo Avendano-Franco'

import pychemia
if pychemia.db.USE_MONGO:
    from pychemia.db import PyChemiaDB
from pychemia import Structure, Composition
import uuid


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

    def __init__(self, name, new=False, composition=None):
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
            self.db = PyChemiaDB(name)
            if new:
                self.db.clean()
        else:
            self.db = {}
        self.members = []
        self.active = []
        self.composition = composition

    @property
    def all_entries(self):
        return [str(entry['_id']) for entry in self.db.entries.find()]

    def get_member_dict(self, imember):
        if pychemia.db.USE_MONGO:
            entry = self.db.entries.find_one({'_id': imember})
        else:
            if imember in self.db:
                entry = self.db[imember]
            else:
                entry = None
        return entry

    def is_evaluated(self, imember):
        entry = self.get_member_dict(imember)
        if entry is not None and entry['properties'] is not None:
            return True
        else:
            return False

    def update_entry(self, imember, structure=None, properties=None, status=None):

        entry = self.get_member_dict(imember)

        if structure is not None:
            entry['structure'] = structure.to_dict()
        if properties is not None:
            entry['properties'] = properties
        if status is not None:
            entry['status'] = status

        if pychemia.db.USE_MONGO:
            self.db.update(imember, entry)
        else:
            self.db[imember] = entry

        return entry

    @staticmethod
    def new_identifier():
        return str(uuid.uuid4())[-12:]

    def new_entry(self, structure, active=True, properties=None):
        entry = {'structure': structure.to_dict(), 'status': {'active': active}, 'properties': properties}

        if pychemia.db.USE_MONGO:
            ident = self.db.entries.insert(entry)
        else:
            ident = self.new_identifier()
            self.db[ident] = entry
        return ident

    def add_random(self):
        """
        Add one random structure to the population
        """
        if self.composition is None:
            raise ValueError('First set a composition')
        structure = Structure.random_cell(self.composition)
        ident = self.new_entry(structure)
        self.actives.append(ident)
        self.members.append(ident)
        return ident

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

    def __repr__(self):
        return str(self.__class__)+'({0.name!r})'.format(self)

    def __str__(self):
        ret = 'Population\n\nName: {0.name!s}\n'.format(self)
        ret += 'Size: %d' % self.db.entries.count()
        return ret

    def __len__(self):
        return self.db.entries.count()

    def set_all_active(self):

        for i in self.all_entries:
            self.active.append(i)

    @property
    def actives(self):
        return self.active

    def add_from_db(self, dbname, sizemax=1):

        comp = Composition(self.composition)
        readdb = PyChemiaDB(dbname)

        index = 0
        for entry in readdb.entries.find({'formula': comp.formula, 'natom': comp.natom}):
            if index < sizemax:
                print 'Adding entry '+str(entry['_id'])+' from '+dbname
                self.new_entry(pychemia.Structure.from_dict(entry))
                index += 1
