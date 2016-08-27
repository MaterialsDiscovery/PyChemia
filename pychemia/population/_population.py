from __future__ import unicode_literals
import json
import numpy as np
from builtins import str
from abc import ABCMeta, abstractmethod
from pychemia import HAS_PYMONGO
from pychemia.utils.computing import deep_unicode

if HAS_PYMONGO:
    from pychemia.db import PyChemiaDB


class Population:
    __metaclass__ = ABCMeta
    """
    General class for all optimization algorithms that uses fixed and blocked
    Generations
    """

    def __init__(self, name, tag, use_mongo=True, direct_evaluation=False):

        name = deep_unicode(name)
        self.tag = tag
        self.pcdb = None
        self.direct_evaluation = direct_evaluation
        if isinstance(name, str):
            self.name = name
            if use_mongo:
                self.pcdb = PyChemiaDB(name)
        else:
            self.name = name.name
            if use_mongo:
                self.pcdb = name

    def __iter__(self):
        return self.pcdb.entries.find()

    def __len__(self):
        return len(self.members)

    def __str__(self):
        ret = ' Population Name:     %s\n' % self.name
        ret += ' Tag:                 %s\n' % self.tag
        ret += ' Members:             %s\n' % len(self)
        return ret

    def disable(self, entry_id):
        self.pcdb.entries.update({'_id': entry_id}, {'$set': {'status.' + self.tag: False}})

    def enable(self, entry_id):
        self.pcdb.entries.update({'_id': entry_id}, {'$set': {'status.' + self.tag: True}})
        if self.direct_evaluation:
            self.evaluate_entry(entry_id)

    def get_values(self, selection):
        ret = {}
        for i in selection:
            ret[i] = self.value(i)
        return ret

    def update_properties(self, entry_id, new_properties):
        self.pcdb.update(entry_id, properties=new_properties)

    def set_in_properties(self, entry_id, field, value):
        return self.pcdb.entries.update_one({'_id': entry_id}, {'$set': {'properties.'+field: value}})

    def get_population_info(self):
        return self.pcdb.db.population_info.find_one({'tag': self.tag})

    def insert_entry(self, entry):
        if 'structure' not in entry:
            entry['structure'] = {}
        if 'properties' not in entry:
            entry['properties'] = {}
        if 'status' not in entry:
            entry['status'] = {}
        return self.pcdb.entries.insert(entry)

    def set_entry(self, entry):
        return self.insert_entry(entry)

    def get_structure(self, entry_id):
        return self.pcdb.get_structure(entry_id)

    def set_structure(self, entry_id, structure):
        return self.pcdb.update(entry_id, structure=structure)

    def unset_properties(self, entry_id):
        return self.pcdb.update(entry_id, properties={})

    def get_entry(self, entry_id, projection=None, with_id=True):
        """
        Return an entry identified by 'entry_id'

        :param with_id:
        :param projection: Insert that projection into the query
        :param entry_id: A database identifier
        :return:
        """
        if not with_id:
            projection['_id'] = 0
        entry = self.pcdb.entries.find_one({'_id': entry_id}, projection)
        return entry

    def ids_sorted(self, selection):
        values = np.array([self.value(i) for i in selection])
        sorted_indices = np.argsort(values)
        return np.array(selection)[sorted_indices]

    def load_json(self, filename):
        filep = open(filename, 'r')
        data = json.load(filep)
        for entry in data:
            self.pcdb.entries.insert(entry)

    def random_population(self, n):
        """
        Create N new random structures to the population

        :param n: (int) The number of new structures
        :return: (list) The identifiers for the new structures
        """
        return [self.add_random() for i in range(n)]

    def replace_failed(self):
        pass

    def save_info(self):
        data = self.pcdb.db.population_info.find_one({'_id': self.tag})
        if data is None:
            data = self.to_dict
            data['_id'] = self.tag
            self.pcdb.db.population_info.insert(data)
        else:
            self.pcdb.db.population_info.update({'_id': self.tag}, self.to_dict)

    def save_json(self, filename):
        ret = []
        for entry_id in self.members:
            ret.append(self.get_entry(entry_id, with_id=False))
        filep = open(filename, 'w')
        json.dump(ret, filep, sort_keys=True, indent=4, separators=(',', ': '))

    def unlock_all(self, name=None):
        for i in self.members:
            self.pcdb.unlock(i, name=name)

    @abstractmethod
    def add_random(self):
        pass

    @abstractmethod
    def check_duplicates(self, ids):
        pass

    @abstractmethod
    def cross(self, ids):
        pass

    @abstractmethod
    def distance(self, entry_id, entry_jd):
        pass

    @abstractmethod
    def evaluate_entry(self, entry_id):
        pass

    @abstractmethod
    def get_duplicates(self, ids):
        pass

    @abstractmethod
    def from_dict(self, population_dict):
        pass

    @abstractmethod
    def is_evaluated(self, entry_id):
        pass

    @abstractmethod
    def move_random(self, entry_id, factor=0.2, in_place=False, kind='move'):
        pass

    @abstractmethod
    def move(self, entry_id, entry_jd, factor=0.2, in_place=False):
        pass

    @abstractmethod
    def new_entry(self, data, active=True):
        pass

    @abstractmethod
    def recover(self):
        pass

    @abstractmethod
    def value(self, entry_id):
        pass

    @abstractmethod
    def str_entry(self, entry_id):
        pass

    @property
    def actives(self):
        return [entry['_id'] for entry in self.pcdb.entries.find({'status.' + self.tag: True}, {'_id': 1})]

    @property
    def actives_evaluated(self):
        return [x for x in self.actives if self.is_evaluated(x)]

    @property
    def actives_no_evaluated(self):
        return [x for x in self.actives if not self.is_evaluated(x)]

    @property
    def evaluated(self):
        return [entry for entry in self.members if self.is_evaluated(entry)]

    @property
    def fraction_evaluated(self):
        ret = np.sum([1 for i in self.actives if self.is_evaluated(i)])
        return float(ret) / len(self.actives)

    @property
    def members(self):
        return [x['_id'] for x in self.pcdb.entries.find({}, {'_id': 1})]

    @property
    def to_dict(self):
        return {'name': self.name, 'tag': self.tag}

    @property
    def best_candidate(self):
        if len(self.evaluated) > 0:
            return self.ids_sorted(self.evaluated)[0]
        else:
            return None

    def refine_progressive(self, entry_id):
        pass
