"""
Includes StructurePopulation
"""

__author__ = 'Guillermo Avendano-Franco'

import uuid
import json
import random
import numpy as np
from fractions import gcd

from pychemia import Composition, Structure, log
from pychemia.db import USE_MONGO

if USE_MONGO:
    from pychemia.db import PyChemiaDB, get_database
from pychemia.analysis import StructureAnalysis, StructureChanger, StructureMatch
from pychemia.utils.mathematics import unit_vector
from pychemia.utils.periodic import atomic_number, covalent_radius


class StructurePopulation():
    def __init__(self, name, composition=None, tag='global', delta=0.3, target_forces=1E-3, value_tol=1E-2,
                 distance_tol=0.3, min_comp_mult=2, max_comp_mult=8, pcdb_source=None):
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
        :return: A new StructurePopulation object
        """
        if composition is not None:
            self.composition = Composition(composition)
        else:
            self.composition = None
        self.delta = delta
        self.tag = tag
        self.target_forces = target_forces
        self.value_tol = value_tol
        self.distance_tol = distance_tol
        self.min_comp_mult = min_comp_mult
        self.max_comp_mult = max_comp_mult
        if isinstance(name, basestring):
            self.name = name
            self.pcdb = PyChemiaDB(name)
        else:
            self.name = name.name
            self.pcdb = name
        self.pcdb_source = pcdb_source
        self.source_blacklist = []

    @property
    def actives(self):
        return [entry['_id'] for entry in self.pcdb.entries.find({'status.' + self.tag: True}, {'_id': 1})]

    @property
    def members(self):
        return [x['_id'] for x in self.pcdb.entries.find({}, {'_id': 1})]

    def __len__(self):
        return len(self.members)

    @property
    def evaluated(self):
        return [entry['_id'] for entry in self.pcdb.entries.find() if self.is_evaluated(entry['_id'])]

    def get_entry(self, entry_id, with_id=True):
        """
        Return an entry identified by 'entry_id'

        :param entry_id: A database identifier
        :return:
        """
        entry = self.pcdb.entries.find_one({'_id': entry_id})
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
        entry_id = self.pcdb.insert(structure=structure, properties=properties, status=status)
        log.debug('Added new entry: %s with tag=%s: %s' % (str(entry_id), self.tag, str(active)))
        return entry_id

    def activate(self, entry_id):
        structure, properties, status = self.pcdb.get_dicts(entry_id)
        status[self.tag] = True
        self.pcdb.update(entry_id, status=status)

    def get_max_force_stress(self, entry_id):
        entry = self.get_entry(entry_id)
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

    def add_random(self, random_probability=0.3):
        """
        Add one random structure to the population
        """
        if self.composition is None:
            raise ValueError('No composition associated to this population')
        comp = self.composition.composition.copy()
        rnd = random.random()
        natom_limit = self.max_comp_mult*self.composition.natom/self.composition.gcd
        condition={'structure.nspecies': self.composition.nspecies,
                   'structure.natom': {'$lte': natom_limit}}

        if self.pcdb_source.entries.find(condition).count() <= len(self.source_blacklist):
            rnd=0

        print 'The random value is ', rnd
        if self.pcdb_source is None or rnd < random_probability or self.composition.nspecies > 1:
            factor = np.random.randint(self.min_comp_mult, self.max_comp_mult + 1)
            for i in comp:
                comp[i] *= factor
            structure = Structure.random_cell(comp)
        else:
            print 'From source...'
            while True:

                entry=None
                condition['properties.spacegroup']=random.randint(1, 230)
                print 'Trying', condition['properties.spacegroup']
                for ientry in self.pcdb_source.entries.find(condition):
                    if ientry['_id'] not in self.source_blacklist:
                        entry=ientry
                        break
                if entry is not None:
                    structure = self.pcdb_source.get_structure(entry['_id'])
                    factor = covalent_radius(self.composition.species[0])/ covalent_radius(structure.species[0])
                    print 'From source: %s Spacegroup: %d Scaling: %7.3f' % (structure.formula,
                                                                             entry['properties']['spacegroup'],
                                                                             factor)
                    structure.set_cell(np.dot(factor*np.eye(3), structure.cell))
                    structure.symbols = structure.natom * self.composition.species
                    self.source_blacklist.append(entry['_id'])
                    break

        return self.new_entry(structure)

    def random_population(self, n):
        """
        Create N new random structures to the population

        :param n: (int) The number of new structures
        :return: (list) The identifiers for the new structures
        """
        return [self.add_random() for i in range(n)]

    def check_duplicates(self, ids):
        ret = {}
        values = np.array([self.value(i) for i in ids])
        if len(values) == 0:
            return ret
        argsort = np.argsort(values)
        diffs = np.ediff1d(values[argsort])

        for i in range(len(values) - 1):
            idiff = diffs[i]
            if idiff < self.value_tol:
                ident1 = ids[argsort[i]]
                ident2 = ids[argsort[i + 1]]
                log.debug('Testing distances between %s and %s' % (str(ident1), str(ident2)))
                distance = self.distance(ident1, ident2)
                # print 'Distance = ', distance
                if distance < self.distance_tol:
                    log.debug('Distance %7.3f < %7.3f' % (distance, self.distance_tol))
                    fx1 = self.value(ident1)
                    fx2 = self.value(ident2)
                    if fx2 < fx1:
                        ret[ident1] = ident2
                    else:
                        ret[ident2] = ident1
        if len(ret) > 0:
            log.debug('Number of duplicates %d' % len(ret))
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

    def distance(self, entry_id, entry_jd, rcut=50):

        fingerprints = {}
        for entry_ijd in [entry_id, entry_jd]:

            if self.pcdb.db.fingerprints.find_one({'_id': entry_ijd}) is None:
                structure = self.get_structure(entry_ijd)
                analysis = StructureAnalysis(structure, radius=rcut)
                x, ys = analysis.fp_oganov()
                fingerprint = {'_id': entry_ijd}
                for k in ys:
                    atomic_number1 = atomic_number(structure.species[k[0]])
                    atomic_number2 = atomic_number(structure.species[k[1]])
                    pair = '%06d' % min(atomic_number1 * 1000 + atomic_number2,
                                        atomic_number2 * 1000 + atomic_number1)
                    fingerprint[pair] = list(ys[k])

                if self.pcdb.db.fingerprints.find_one({'_id': entry_ijd}) is None:
                    self.pcdb.db.fingerprints.insert(fingerprint)
                else:
                    self.pcdb.db.fingerprints.update({'_id': entry_ijd}, fingerprint)
                fingerprints[entry_ijd] = fingerprint
            else:
                fingerprints[entry_ijd] = self.pcdb.db.fingerprints.find_one({'_id': entry_ijd})

        dij = []

        for pair in fingerprints[entry_id]:
            if pair in fingerprints[entry_jd] and pair != '_id':
                uvect1 = unit_vector(fingerprints[entry_id][pair])
                uvect2 = unit_vector(fingerprints[entry_jd][pair])
                dij.append(0.5 * (1.0 - np.dot(uvect1, uvect2)))
        return float(np.mean(dij))

    def add_from_db(self, db_settings, sizemax=1):
        if self.composition is None:
            raise ValueError('No composition associated to this population')
        comp = Composition(self.composition)
        readdb = get_database(db_settings)

        index = 0
        for entry in readdb.entries.find({'structure.formula': comp.formula,
                                          'structure.natom': {'$lte': self.min_comp_mult * comp.natom,
                                                              '$gte': self.max_comp_mult * comp.natom}}):
            if index < sizemax:
                print 'Adding entry ' + str(entry['_id']) + ' from ' + readdb.name
                self.new_entry(readdb.get_structure(entry['_id']))
                index += 1

    def add_modified(self, entry_id):

        structure = self.pcdb.get_structure(entry_id)
        changer = StructureChanger(structure)
        changer.random_change(self.delta)
        new_structure = changer.new_structure
        return self.new_entry(new_structure)

    def disable(self, entry_id):
        self.pcdb.entries.update({'_id': entry_id}, {'$set': {'status.'+self.tag: False}})

    def enable(self, entry_id):
        self.pcdb.entries.update({'_id': entry_id}, {'$set': {'status.'+self.tag: True}})

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
            self.pcdb.entries.insert(entry)

    def move_random(self, entry_id, factor=0.2, in_place=False):
        structure = self.get_structure(entry_id)

        changer = StructureChanger(structure=structure)
        changer.random_move_many_atoms(epsilon=factor)
        if in_place:
            return self.pcdb.update(entry_id, structure=changer.new_structure)
        else:
            return self.new_entry(changer.new_structure, active=False)

    def move(self, entry_id, entry_jd, factor=0.2, in_place=False):
        """
        Moves entry_id in the direction of entry_jd
        If in_place is True the movement occurs on the
        same address as entry_id

        :param entry_id:
        :param entry_jd:
        :param in_place:
        :return:
        """
        structure_mobile = self.get_structure(entry_id)
        structure_target = self.get_structure(entry_jd)

        if structure_mobile.natom != structure_target.natom:
            # Moving structures with different number of atoms is only implemented for smaller structures moving
            # towards bigger ones by making a super-cell and only if their size is smaller that 'max_comp_mult'

            mult1 = structure_mobile.get_composition().gcd
            mult2 = structure_target.get_composition().gcd
            lcd = mult1 * mult2 / gcd(mult1, mult2)
            if lcd > self.max_comp_mult:
                # The resulting structure is bigger than the limit
                # cannot move
                if not in_place:
                    return self.new_entry(structure_mobile)
                else:
                    return entry_id

        # We will move structure1 in the direction of structure2
        match = StructureMatch(structure_target, structure_mobile)
        match.match_atoms()
        displacements = match.reduced_displacement()

        new_reduced = match.structure2.reduced + factor * displacements
        new_cell = match.structure2.cell
        new_symbols = match.structure2.symbols
        new_structure = Structure(reduced=new_reduced, symbols=new_symbols, cell=new_cell)
        if in_place:
            return self.pcdb.update(entry_id, structure=new_structure)
        else:
            return self.new_entry(new_structure)

    def __str__(self):
        ret = ' Structure Population\n\n'
        ret += ' Name:               %s\n' % self.name
        ret += ' Tag:                %s\n' % self.tag
        ret += ' Delta:              %7.2E\n' % self.delta
        ret += ' Target-Forces:      %7.2E\n' % self.target_forces
        ret += ' Value tolerance:    %7.2E\n' % self.value_tol
        ret += ' Distance tolerance: %7.2E\n\n' % self.distance_tol
        if self.composition is not None:
            ret += ' Composition:                  %s\n' % self.composition.formula
            ret += ' Minimal composition multiplier: %d\n' % self.min_comp_mult
            ret += ' Maximal composition multiplier: %d\n\n' % self.max_comp_mult
        else:
            ret += '\n'
        ret += ' Members:            %d\n' % len(self.members)
        ret += ' Actives:            %d\n' % len(self.actives)
        ret += ' Evaluated:          %d\n' % len(self.evaluated)
        return ret

    def unlock_all(self, name=None):
        for i in self.members:
            self.pcdb.unlock(i, name=name)

    def ids_sorted(self, selection):
        values = np.array([self.value(i) for i in selection])
        sorted_indices = np.argsort(values)
        return np.array(selection)[sorted_indices]

    def get_values(self, selection):
        ret = {}
        for i in selection:
            ret[i] = self.value(i)
        return ret

    def value(self, entry_id):
        entry = self.get_entry(entry_id)
        structure = self.get_structure(entry_id)
        if 'properties' not in entry:
            log.debug('This entry has no properties %s' % str(entry['_id']))
            return None
        elif entry['properties'] is None:
            return None
        elif 'energy' not in entry['properties']:
            log.debug('This entry has no energy in properties %s' % str(entry['_id']))
            return None
        else:
            return entry['properties']['energy'] / structure.get_composition().gcd

    def to_dict(self):
        return {'name': self.name,
                'tag': self.tag,
                'delta': self.delta,
                'target_forces': self.target_forces,
                'value_tol': self.value_tol,
                'distance_tol': self.distance_tol}

    @staticmethod
    def from_dict(population_dict):
        return StructurePopulation(name=population_dict['name'],
                                   tag=population_dict['tag'],
                                   delta=population_dict['delta'],
                                   target_forces=population_dict['target_forces'],
                                   value_tol=population_dict['value_tol'],
                                   distance_tol=population_dict['distance_tol'])

    def save_info(self):
        self.pcdb.db.population_info.insert(self.to_dict())

    def __iter__(self):
        return self.pcdb.entries.find()
