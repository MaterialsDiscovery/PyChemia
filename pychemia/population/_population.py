"""
Includes StructurePopulation
"""

__author__ = 'Guillermo Avendano-Franco'

import uuid
import json
import numpy as np

from pychemia import Composition, Structure, log
from pychemia.db import USE_MONGO
if USE_MONGO:
    from pychemia.db import PyChemiaDB
from pychemia.analysis import StructureAnalysis, StructureChanger
from pychemia.utils.mathematics import unit_vector


class StructurePopulation():
    def __init__(self, name, composition, tag='global', delta=0.1, new=False, target_forces=1E-3):
        """
        Defines a population of PyChemia Structures,

        The 'name' of the database is used to create the MongoDB database and the structures are
        uniform in composition. A specific 'tag' could be attached to differentiate
        the other instances running concurrently. The 'delta' argument is the scaling
        factor for changers and mixers. In the case of populations supported on
        PyChemia databases the 'new' will erase the database

        :param name: The name of the population. ie the name of the database
        :param composition: The composition uniform for all the members
        :param tag: A tag to differentiate different instances running concurrently
        :param delta: The parameter to scale the changers and mixers
        :param new: If true the database will be erased
        :return: A new StructurePopulation object
        """
        self.composition = Composition(composition)
        self.delta = delta
        self.name = name
        self.tag = tag
        self.target_forces = target_forces

        self._members = []
        self._actives = []
        self._evaluated = []
        if USE_MONGO:
            self.db = PyChemiaDB(name)
            if new:
                self.db.clean()
            else:
                self.check_tags()
        else:
            self.db = {}

    def check_tags(self):
        if USE_MONGO:
            for entry in self.db.entries.find():
                if self.tag not in entry['status']:
                    new_status = entry['status'].copy()
                    new_status[self.tag] = False
                    self.update_entry(entry['_id'], status=new_status)
        else:
            for entry in self.db:
                if self.tag not in self.db[entry]['status']:
                    new_status = db[entry]['status'].copy()
                    new_status[self.tag] = False
                    self.update_entry(entry, status=new_status)

    @property
    def actives(self):
        self.check_tags()
        if USE_MONGO:
            self._actives = [entry['_id'] for entry in self.db.entries.find() if
                             self.tag in entry['status'] and entry['status'][self.tag]]
        return self._actives

    @property
    def members(self):
        if USE_MONGO:
            self._members = [entry['_id'] for entry in self.db.entries.find()]
        return self._members

    @property
    def evaluated(self):
        if USE_MONGO:
            self._evaluated = [entry['_id'] for entry in self.db.entries.find() if self.is_evaluated(entry['_id'])]
        return self._evaluated

    def get_member_dict(self, imember, with_id=True):
        """
        Return an entry identified by 'imember'

        :param imember: A database identifier
        :return:
        """
        if USE_MONGO:
            entry = self.db.entries.find_one({'_id': imember})
            if not with_id:
                entry.pop('_id')
        else:
            if imember in self.db:
                entry = self.db[imember]
                if with_id:
                    entry['_id'] = imember
            else:
                entry = None
        return entry

    def get_structure(self, imember):
        entry = self.get_member_dict(imember)
        return Structure.from_dict(entry['structure'])

    def get_properties(self, imember):
        entry = self.get_member_dict(imember)
        return Structure.from_dict(entry['properties'])

    def get_status(self, imember):
        entry = self.get_member_dict(imember)
        return Structure.from_dict(entry['status'])

    def update_entry(self, imember, structure=None, properties=None, status=None):

        entry = self.get_member_dict(imember)

        if structure is not None:
            entry['structure'] = structure.to_dict()
        if properties is not None:
            entry['properties'] = properties
        if status is not None:
            entry['status'] = status

        if USE_MONGO:
            self.db.update(imember, entry)
        else:
            self.db[imember] = entry
        return imember

    @staticmethod
    def new_identifier():
        return str(uuid.uuid4())[-12:]

    def new_entry(self, structure, active=True, properties=None):
        entry = {'structure': structure.to_dict(), 'status': {self.tag: active}, 'properties': properties}

        if USE_MONGO:
            ident = self.db.entries.insert(entry)
        else:
            ident = self.new_identifier()
            self.db[ident] = entry
            self._actives.append(ident)
            self._members.append(ident)

        return ident

    def is_evaluated(self, imember):
        entry = self.get_member_dict(imember)
        if entry is not None and entry['properties'] is not None:
            if 'energy' not in entry['properties']:
                return False
            if 'stress' not in entry['properties']:
                return False
            if 'forces' not in entry['properties']:
                return False
            if 'forces' in entry['properties'] and np.max(np.abs(np.array(entry['properties']['forces'], dtype=float).flatten())) > self.target_forces:
                return False
            if 'stress' in entry['properties'] and np.max(np.abs(np.array(entry['properties']['stress'], dtype=float).flatten())) > self.target_forces:
                return False
            return True
        else:
            return False

    def add_random(self):
        """
        Add one random structure to the population
        """
        if self.composition is None:
            raise ValueError('First set a composition')
        structure = Structure.random_cell(self.composition)
        ident = self.new_entry(structure)
        return ident

    def random_population(self, n):
        """
        Create N new random structures to the population

        :param n: (int) The number of new structures
        :return: (list) The identifiers for the new structures
        """
        ret = []
        for i in range(n):
            ret.append(self.add_random())
        return ret

    def value(self, imember):
        entry = self.get_member_dict(imember)
        return entry['properties']['energy']

    def check_duplicates(self, value_tol=1E-2, distance_tol=0.3):
        ret = []
        ids = self.actives
        values = np.array([self.value(i) for i in self.actives if i in self.evaluated])
        print 'Values= ', sorted(values)
        if len(values) == 0:
            return ret
        # else:
        #    print 'Values', values
        argsort = np.argsort(values)
        diffs = np.ediff1d(values[argsort])
        #print 'Values     ', values[argsort]
        #print 'Differences', diffs
        for i in range(len(values) - 1):
            idiff = diffs[i]
            if idiff < value_tol:
                ident1 = ids[argsort[i]]
                ident2 = ids[argsort[i + 1]]
                print 'Testing distances between ', ident1, ' and ', ident2
                distance = self.distance(ident1, ident2)
                print 'Distance = ', distance
                if distance < distance_tol:
                    print "Very small distance ", distance
                    fx1 = self.value(ident1)
                    fx2 = self.value(ident2)
                    if fx2 < fx1:
                        ret.append(ident1)
                    else:
                        ret.append(ident2)
        if len(ret) > 0:
            print 'Duplicates', ret
        else:
            print 'No duplicates'
        return ret

    def distance_matrix(self):

        members = self.members
        analysis = {}
        ret = np.zeros((len(members), len(members)))
        for i in members:
            analysis[i] = StructureAnalysis(self.get_structure(i))

        for i in range(len(members)):
            log.debug('Fingerprint for %s' % str(members[i]))
            x1, y1_dict = analysis[members[i]].fp_oganov()
            for j in range(i, len(members)):
                log.debug('Fingerprint for %s' % str(members[j]))
                x2, y2_dict = analysis[members[j]].fp_oganov()

                dij = []
                for x in y1_dict:
                    uvect1 = unit_vector(y1_dict[x])
                    uvect2 = unit_vector(y2_dict[x])
                    dij.append(0.5 * (1.0 - np.dot(uvect1, uvect2)))
                ret[i, j] = np.mean(dij)
                ret[j, i] = ret[i, j]
        return ret

    def diff_values_matrix(self):

        members = self.members
        ret = np.zeros((len(members), len(members)))

        for i in range(len(members)):
            for j in range(i, len(members)):

                if self.value(members[i]) is not None and self.value(members[j]) is not None:
                    ret[i, j] = np.abs(self.value(members[i])-self.value(members[j]))
                else:
                    ret[i, j] = float('nan')
                ret[j, i] = ret[i, j]
        return ret

    def distance(self, imember, jmember, rcut=50):
        struct1 = Structure.from_dict(self.get_member_dict(imember)['structure'])
        struct2 = Structure.from_dict(self.get_member_dict(jmember)['structure'])
        analysis1 = StructureAnalysis(struct1, radius=rcut)
        analysis2 = StructureAnalysis(struct2, radius=rcut)
        x1, y1_dict = analysis1.fp_oganov()
        x2, y2_dict = analysis2.fp_oganov()
        # print len(x1)
        assert (len(x1) == len(x2))
        #print np.dot(unit_vector(x1), unit_vector(x2))
        assert (len(y1_dict) == len(y2_dict))
        dij = []
        for i in y1_dict:
            uvect1 = unit_vector(y1_dict[i])
            uvect2 = unit_vector(y2_dict[i])
            dij.append(0.5 * (1.0 - np.dot(uvect1, uvect2)))
        return np.mean(dij)

    def add_from_db(self, dbname, sizemax=1):

        comp = Composition(self.composition)
        readdb = PyChemiaDB(dbname)

        index = 0
        for entry in readdb.entries.find({'formula': comp.formula, 'natom': comp.natom}):
            if index < sizemax:
                print 'Adding entry ' + str(entry['_id']) + ' from ' + dbname
                self.new_entry(Structure.from_dict(entry))
                index += 1

    def add_modified(self, ident):

        structure = self.get_structure(ident)

        changer = StructureChanger(structure)
        changer.random_change(self.delta)

        new_structure = changer.new_structure

        new_ident = self.new_entry(new_structure)
        return new_ident

    def disable(self, ident):
        if ident not in self.actives:
            raise ValueError(str(ident) + ' not in actives')
        entry = self.get_member_dict(ident)
        entry['status'][self.tag] = False
        structure = Structure.from_dict(entry['structure'])
        self.update_entry(ident, structure=structure, properties=entry['properties'], status=entry['status'])
        if not USE_MONGO:
            self._actives.remove(ident)

    @property
    def fraction_evaluated(self):
        print 'Evaluated: ', self.evaluated
        print 'Actives:', self.actives
        ret = np.sum([1 for i in self.actives if i in self.evaluated])
        return float(ret) / len(self.actives)

    @property
    def actives_no_evaluated(self):
        ret = []
        for i in self.actives:
            if i not in self.evaluated:
                ret.append(i)
        return ret

    def save(self):
        ret = []
        for imember in self.members:
            ret.append(self.get_member_dict(imember, with_id=False))
        filep = open(self.name + '.pop', 'w')
        json.dump(ret, filep, sort_keys=True, indent=4, separators=(',', ': '))

    def member_str(self, imember):
        data = self.get_member_dict(imember)
        ret = imember
        return ret

    def move(self, imember, jmember, in_place=False):
        """
        Moves imember in the direction of jmember
        If in_place is True the movement occurs on the
        same address as imember

        :param imember:
        :param jmember:
        :param in_place:
        :return:
        """
        structure1 = self.get_structure(imember)
        structure2 = self.get_structure(jmember)

        assert (structure1.symbols == structure2.symbols)

        new_reduced = np.zeros((structure1.natom, 3))
        for i in range(structure1.natom):
            x1 = structure1.reduced[i]
            x2 = structure2.reduced[i]
            uvector = (x2 - x1) / np.linalg.norm(x2 - x1)
            new_reduced[i] = structure1.reduced[i] + self.delta * uvector
        new_cell = 0.5 * (structure1.cell + structure2.cell)
        new_structure = Structure(reduced=new_reduced, symbols=structure1.symbols, cell=new_cell)
        if not in_place:
            new_ident = self.new_entry(new_structure)
        else:
            new_ident = self.update_entry(imember, new_structure)
        return new_ident

    def complete_actives(self, n):
        difference = n - len(self.actives)
        if difference > 0:
            print 'Adding ', difference, 'new random structures'
            self.random_population(difference)

    def __str__(self):
        ret = ' Structure Population\n\n'
        ret += ' Name:        ' + self.name + '\n'
        ret += ' Composition: ' + str(self.composition) + '\n'
        ret += ' Delta:       ' + str(self.delta) + '\n'
        ret += ' Tag:         ' + self.tag + '\n'
        ret += '\n'
        ret += ' Members:         ' + str(len(self.members)) + '\n'
        ret += ' Actives:         ' + str(len(self.actives)) + '\n'
        ret += ' Evaluated:       ' + str(len(self.evaluated)) + '\n'
        return ret

class ObjectiveFunction():
    def __init__(self):
        self.population = None

    def initialize(self, population):
        self.population = population

    def ids_sorted(self, selection):
        values = np.array([self.population.value(i) for i in selection])
        argsort = np.argsort(values)
        return np.array(selection)[argsort]

    def get_values(self, selection):
        ret = {}
        for i in selection:
            ret[i] = self.population.value(i)
        return ret
