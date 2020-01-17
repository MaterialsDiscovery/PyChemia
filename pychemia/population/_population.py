
import json
import numpy as np
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

    def __init__(self, name, tag, use_mongo=True, direct_evaluation=False, distance_tolerance=0.1):

        name = deep_unicode(name)
        self.tag = tag
        self.pcdb = None
        self.direct_evaluation = direct_evaluation
        self.distance_tolerance = distance_tolerance
        if isinstance(name, str):
            self.name = name
            if use_mongo:
                self.pcdb = PyChemiaDB(name)
        else:
            self.name = name.name
            if use_mongo:
                self.pcdb = name

    def __iter__(self):
        if self.tag != 'global':
            return self.pcdb.entries.find({'status.tag': self.tag})
        else:
            return self.pcdb.entries.find({})

    def __len__(self):
        return len(self.members)

    def __str__(self):
        ret = '[%s] Database: %s\n' % (self.tag, self.name)
        ret += '[%s] Tag:      %s\n' % (self.tag, self.tag)
        ret += '[%s] Members:  %s\n' % (self.tag, len(self))
        return ret

    def disable(self, entry_id):
        self.pcdb.entries.update_one({'_id': entry_id}, {'$set': {'status.' + self.tag: False}})

    def enable(self, entry_id):
        self.pcdb.entries.update_one({'_id': entry_id}, {'$set': {'status.' + self.tag: True}})
        if self.direct_evaluation:
            self.evaluate_entry(entry_id)

    def get_values(self, selection):
        ret = {}
        for i in selection:
            ret[i] = self.value(i)
        return ret

    def clean(self):
        self.pcdb.clean()

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
        result = self.pcdb.entries.insert_one(entry)
        return result.inserted_id

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
            self.pcdb.entries.insert_one(entry)

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
            self.pcdb.db.population_info.insert_one(data)
        else:
            self.pcdb.db.population_info.replace_one({'_id': self.tag}, self.to_dict)

    def save_json(self, filename):
        ret = []
        for entry_id in self.members:
            ret.append(self.get_entry(entry_id, with_id=False))
        filep = open(filename, 'w')
        json.dump(ret, filep, sort_keys=True, indent=4, separators=(',', ': '))

    def unlock_all(self, name=None):
        for i in self.members:
            self.pcdb.unlock(i, name=name)

    def distance_matrix(self, ids):

        ret = np.zeros((len(ids), len(ids)))

        for i in range(len(ids) - 1):
            for j in range(i, len(ids)):
                ret[i, j] = self.distance(ids[i], ids[j])
                ret[j, i] = ret[i, j]
        return ret

    def get_duplicates(self, ids, tolerance=None, fast=True):
        """
        For a given list of identifiers 'ids' checks the values for the function 'distance' and return a dictionary
          where each key is the identifier of a duplicate candidate and the value is a list of identifiers considered
          equivalents to it.


        :param tolerance: When the value of "distance" is lower than the tolerance the candidates are considered
                            equivalent
        :param fast: (bool) Once one structure is considered duplicated do not perform a search for more equivalent
                            candidates, it means that the values of the dictionary are lists with one one member.
                            Otherwise a extended search is done resulting on a complete N^2 search for all candidates
                            versus each other.
        :return:
        :param ids:  List of identifiers for wich the check will be performed
        :return:
        """
        if tolerance is None:
            tolerance = self.distance_tolerance
        ids = self.ids_sorted(ids)
        ret = {}
        for i in range(len(ids)-1):
            if fast and ids[i] in ret:
                continue
            for j in range(i+1, len(ids)):
                if self.distance(ids[i], ids[j]) < tolerance:
                    if fast:
                        ret[ids[j]] = ids[i]
                    else:
                        if ids[j] in ret:
                            ret[ids[j]].append(ids[i])
                        else:
                            ret[ids[j]] = [ids[i]]
        return ret

    @abstractmethod
    def add_random(self):
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
        if len(self.actives) != 0:
            return float(ret) / len(self.actives)
        else:
            return 0

    @property
    def members(self):
        if self.tag != 'global':
            return [x['_id'] for x in self.pcdb.entries.find({'status.tag': self.tag}, {'_id': 1})]
        else:
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
