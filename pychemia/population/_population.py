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
    from pychemia.db import PyChemiaDB, get_database
from pychemia.analysis import StructureAnalysis, StructureChanger
from pychemia.utils.mathematics import unit_vector


class StructurePopulation():
    def __init__(self, name, composition, tag='global', delta=0.1, target_forces=1E-3, value_tol=1E-2,
                 distance_tol=0.3, min_comp_mult=2, max_comp_mult=8):
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
        self.value_tol = value_tol
        self.distance_tol = distance_tol
        self.min_comp_mult = min_comp_mult
        self.max_comp_mult = max_comp_mult
        self.db = PyChemiaDB(name)

    @property
    def actives(self):
        return [entry['_id'] for entry in self.db.entries.find({'status.' + self.tag: True}, {'_id': 1})]

    @property
    def members(self):
        return [x['_id'] for x in self.db.entries.find({}, {'_id': 1})]

    @property
    def evaluated(self):
        return [entry['_id'] for entry in self.db.entries.find() if self.is_evaluated(entry['_id'])]

    def get_entry(self, entry_id, with_id=True):
        """
        Return an entry identified by 'entry_id'

        :param entry_id: A database identifier
        :return:
        """
        entry = self.db.entries.find_one({'_id': entry_id})
        if entry is not None and not with_id:
            entry.pop('_id')
        return entry

    def get_structure(self, entry_id):
        entry = self.get_entry(entry_id)
        return Structure.from_dict(entry['structure'])

    @staticmethod
    def new_identifier():
        return str(uuid.uuid4())[-12:]

    def new_entry(self, structure, active=True):
        properties = {'forces': None, 'stress': None, 'energy': None}
        status = {self.tag: active}
        return self.db.insert(structure=structure, properties=properties, status=status)

    def get_max_force_stress(self, imember):
        entry = self.get_entry(imember)
        if entry is not None and entry['properties'] is not None:
            properties = entry['properties']
            if 'forces' not in properties or 'stress' not in properties:
                forces = None
                stress = None
            elif properties['forces'] is None or properties['stress'] is None:
                forces = None
                stress = None
            else:
                forces = np.max(np.abs(np.array(entry['properties']['forces'], dtype=float).flatten()))
                stress = np.max(np.abs(np.array(entry['properties']['stress'], dtype=float).flatten()))
        else:
            forces = None
            stress = None
        return forces, stress

    def is_evaluated(self, entry_id):
        max_force, max_stress = self.get_max_force_stress(entry_id)
        if max_force is None or max_stress is None:
            return False
        elif max_force < self.target_forces and max_stress < self.target_forces:
            return True
        else:
            return False

    def add_random(self):
        """
        Add one random structure to the population
        """
        factor = np.random.randint(self.min_comp_mult, self.max_comp_mult + 1)
        comp = self.composition.composition.copy()
        for i in comp:
            comp[i] *= factor

        structure = Structure.random_cell(comp)
        return self.new_entry(structure)

    def random_population(self, n):
        """
        Create N new random structures to the population

        :param n: (int) The number of new structures
        :return: (list) The identifiers for the new structures
        """
        return [self.add_random() for i in range(n)]

    def check_duplicates(self):
        ret = []
        ids = self.actives_evaluated
        values = np.array([self.value(i) for i in ids])
        log.debug('Values: %s' % str(sorted(values)))
        if len(values) == 0:
            return ret
        argsort = np.argsort(values)
        diffs = np.ediff1d(values[argsort])

        log.debug('Differences: %s' % str(diffs))
        for i in range(len(values) - 1):
            idiff = diffs[i]
            if idiff < self.value_tol:
                ident1 = ids[argsort[i]]
                ident2 = ids[argsort[i + 1]]
                log.debug('Testing distances between %s and %s' % (str(ident1), str(ident2)))
                distance = self.distance(ident1, ident2)
                print 'Distance = ', distance
                if distance < self.distance_tol:
                    log.debug('Distance %7.3f < %7.3f' % (distance, self.distance_tol))
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

    def distance_matrix(self, radius=20):

        members = self.members
        analysis = {}
        ret = np.zeros((len(members), len(members)))
        for i in members:
            analysis[i] = StructureAnalysis(self.get_structure(i), radius=radius)

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
                    ret[i, j] = np.abs(self.value(members[i]) - self.value(members[j]))
                else:
                    ret[i, j] = float('nan')
                ret[j, i] = ret[i, j]
        return ret

    def distance(self, imember, jmember, rcut=50):
        struct1 = self.get_structure(imember)
        struct2 = self.get_structure(jmember)
        analysis1 = StructureAnalysis(struct1, radius=rcut)
        analysis2 = StructureAnalysis(struct2, radius=rcut)
        x1, y1_dict = analysis1.fp_oganov()
        x2, y2_dict = analysis2.fp_oganov()
        assert (len(x1) == len(x2))
        assert (len(y1_dict) == len(y2_dict))
        dij = []
        for i in y1_dict:
            uvect1 = unit_vector(y1_dict[i])
            uvect2 = unit_vector(y2_dict[i])
            dij.append(0.5 * (1.0 - np.dot(uvect1, uvect2)))
        return float(np.mean(dij))

    def add_from_db(self, db_settings, sizemax=1):

        comp = Composition(self.composition)
        readdb = get_database(db_settings)

        index = 0
        for entry in readdb.entries.find({'structure.formula': comp.formula,
                                          'structure.natom': {'$lte': self.min_comp_mult * comp.natom,
                                                              '$gte': self.max_comp_mult * comp.natom}}):
            if index < sizemax:
                print 'Adding entry ' + str(entry['_id']) + ' from ' + dbname
                self.new_entry(readdb.get_structure(entry['_id']))
                index += 1

    def add_modified(self, entry_id):

        structure = self.db.get_structure(entry_id)
        changer = StructureChanger(structure)
        changer.random_change(self.delta)
        new_structure = changer.new_structure
        return self.new_entry(new_structure)

    def disable(self, entry_id):
        entry = self.get_entry(entry_id)
        status = None
        if 'status' in entry:
            status = entry['status']
        if status is None:
            status = {}
        status[self.tag] = False
        self.db.update(entry_id, status=status)

    def enable(self, entry_id):
        entry = self.get_entry(entry_id)
        status = None
        if 'status' in entry:
            status = entry['status']
        if status is None:
            status = {}
        status[self.tag] = True
        self.db.update(entry_id, status=status)

    @property
    def fraction_evaluated(self):
        ret = np.sum([1 for i in self.actives if self.is_evaluated(i)])
        return float(ret) / len(self.actives)

    @property
    def actives_no_evaluated(self):
        return [x for x in self.actives if not self.is_evaluated(x)]

    @property
    def actives_evaluated(self):
        return [x for x in self.actives if self.is_evaluated(x)]

    def save_json(self, filename):
        ret = []
        for entry_id in self.members:
            ret.append(self.get_entry(entry_id, with_id=False))
        filep = open(filename, 'w')
        json.dump(ret, filep, sort_keys=True, indent=4, separators=(',', ': '))

    def load_json(self, filename):
        filep = open(filename, 'r')
        data = json.load(filep)
        for entry in data:
            self.db.entries.insert(entry)

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
            return self.new_entry(new_structure)
        else:
            return self.db.update(imember, structure=new_structure)

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

    def unlock_all(self, name=None):
        for i in self.members:
            self.db.unlock(i, name=name)

    def ids_sorted(self, selection):
        values = np.array([self.value(i) for i in selection])
        argsort = np.argsort(values)
        return np.array(selection)[argsort]

    def get_values(self, selection):
        ret = {}
        for i in selection:
            ret[i] = self.value(i)
        return ret

    def value(self, imember):
        entry = self.get_entry(imember)
        struct = self.get_structure(imember)
        if 'properties' not in entry:
            log.debug('This entry has no properties %s' % str(entry['_id']))
            return None
        elif entry['properties'] is None:
            return None
        elif 'energy' not in entry['properties']:
            log.debug('This entry has no energy in properties %s' % str(entry['_id']))
            return None
        else:
            return entry['properties']['energy'] / struct.get_composition().gcd
