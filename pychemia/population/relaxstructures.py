import random
import uuid
from math import gcd
import numpy as np
from ._population import Population
from pychemia import Composition, Structure, pcm_log
from pychemia.analysis import StructureAnalysis, StructureChanger, StructureMatch
from pychemia.analysis.splitting import SplitMatch
from pychemia.utils.mathematics import unit_vector
from pychemia.utils.periodic import atomic_number, covalent_radius
from pymongo import ASCENDING
from pychemia.db import get_database
from pychemia.crystal import CrystalSymmetry


class RelaxStructures(Population):

    def evaluate_entry(self, entry_id):
        pass

    def __init__(self, name, composition=None, tag='global', target_forces=1E-3, value_tol=1E-2,
                 distance_tolerance=0.3, min_comp_mult=2, max_comp_mult=8, pcdb_source=None, pressure=0.0,
                 target_stress=None, target_diag_stress=None, target_nondiag_stress=None):
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
        :return: A new StructurePopulation object
        """
        if composition is not None:
            self.composition = Composition(composition)
        else:
            self.composition = None
        self.tag = tag
        self.target_forces = target_forces
        self.value_tol = value_tol
        self.min_comp_mult = min_comp_mult
        self.max_comp_mult = max_comp_mult
        self.pcdb_source = pcdb_source
        self.pressure = pressure
        if target_stress is None:
            self.target_stress = target_forces
        else:
            self.target_stress = target_stress
        if target_diag_stress is None:
            self.target_diag_stress = self.target_stress
        else:
            self.target_diag_stress = target_diag_stress
        if target_diag_stress is None:
            self.target_nondiag_stress = self.target_stress
        else:
            self.target_nondiag_stress = target_nondiag_stress
        self.name = name
        Population.__init__(self, name, tag, distance_tolerance=distance_tolerance)

        if self.pcdb_source is not None:
            self.sources = {}
            for i in range(min_comp_mult, max_comp_mult+1):
                self.sources[i] = []
                for entry in self.pcdb_source.entries.find({'structure.natom': i*self.composition.natom,
                                                            'structure.nspecies': self.composition.nspecies},
                                                           {'_id': 1}):
                    self.sources[i].append(entry['_id'])

    def recover(self):
        data = self.get_population_info()
        if data is not None:
            self.distance_tolerance = data['distance_tol']
            self.value_tol = data['value_tol']
            self.name = data['name']
            self.target_forces = data['target_forces']

    def get_structure(self, entry_id):
        entry = self.get_entry(entry_id)
        return Structure.from_dict(entry['structure'])

    @staticmethod
    def new_identifier():
        return str(uuid.uuid4())[-12:]

    def new_entry(self, structure, active=True):
        properties = {'forces': None, 'stress': None, 'energy': None}
        status = {self.tag: active, 'tag': self.tag}
        entry = {'structure': structure.to_dict, 'properties': properties, 'status': status}
        entry_id = self.insert_entry(entry)
        pcm_log.debug('Added new entry: %s with tag=%s: %s' % (str(entry_id), self.tag, str(active)))
        return entry_id

    def get_max_force_stress(self, entry_id):
        entry = self.get_entry(entry_id, projection={'properties': 1})
        max_force = None
        max_diag_stress = None
        max_nondiag_stress = None
        if entry is not None and entry['properties'] is not None:
            properties = entry['properties']
            if 'forces' in properties and 'stress' in properties:
                if properties['forces'] is not None and properties['stress'] is not None:
                    forces = np.array(entry['properties']['forces'])
                    stress = np.array(entry['properties']['stress'])
                    max_force = np.max(np.apply_along_axis(np.linalg.norm, 1, forces))
                    max_diag_stress = np.max(np.abs(stress[:3]))
                    max_nondiag_stress = np.max(np.abs(stress[4:]))

        return max_force, max_diag_stress, max_nondiag_stress

    def is_evaluated(self, entry_id):
        max_force, max_diag_stress, max_nondiag_stress = self.get_max_force_stress(entry_id)
        if max_force is None or max_diag_stress is None or max_nondiag_stress is None:
            return False
        elif max_force < self.target_forces and max_diag_stress < self.target_diag_stress + self.pressure:
            if max_nondiag_stress < self.target_nondiag_stress:
                return True
            else:
                return False
        else:
            return False

    def add_random(self, random_probability=0.3):
        """
        Add one random structure to the population
        """
        entry_id = None
        structure = Structure()
        if self.composition is None:
            raise ValueError('No composition associated to this population')
        factor = np.random.randint(self.min_comp_mult, self.max_comp_mult + 1)
        comp = self.composition.composition.copy()
        # print("Initial composition: %s" % comp)
        # print(Composition(comp))
        # print(Composition(comp).symbols)
        for i in comp:
            comp[i] *= factor
        new_comp = Composition(comp)

        while True:
            rnd = random.random()
            condition = {'structure.nspecies': new_comp.nspecies,
                         'structure.natom': new_comp.natom}
            if self.pcdb_source is None:
                rnd = 0
            elif len(self.sources[factor]) == 0:
                rnd = 0
            if self.pcdb_source is None or rnd < random_probability:
                pcm_log.debug('Random Structure')
                structure = Structure.random_cell(new_comp, method='stretching', stabilization_number=5, nparal=5,
                                                  periodic=True)
                break
            else:
                pcm_log.debug('From source')
                entry_id = self.sources[factor][np.random.randint(0, len(self.sources[factor]))]
                structure = self.pcdb_source.get_structure(entry_id)
                print("chosen structure from database =", structure)
                sym = CrystalSymmetry(structure)

                scale_factor = float(np.max(covalent_radius(new_comp.species)) /
                                     np.max(covalent_radius(structure.species)))
                reduce_scale = scale_factor ** (1. / 3)    # WIH
                msg = 'Mult: %d natom: %d From source: %s Spacegroup: %d Scaling: %7.3f'
                print(msg % (factor, structure.natom, structure.formula, sym.number(), scale_factor))
                # structure.set_cell(np.dot(scale_factor * np.eye(3), structure.cell)) # WIH
                structure.set_cell(np.dot(reduce_scale * np.eye(3), structure.cell))  # WIH
                print("symbols before change = ", structure.symbols)
                structure.symbols = new_comp.symbols
                print("symbols after change = ", structure.symbols)
                self.sources[factor].remove(entry_id)
                break

        return self.new_entry(structure), entry_id

    def check_duplicates(self, ids):
        """
        Computes duplicate structures measuring its distance when their value is larger than value_tol.
        If the distance is lower than 'distance_tol' the structures will be cosidered as duplicates.

        :param ids:
        :return: (dict) Dictionary of duplicates, the keys are the ids of the duplicates and the value is the structure
                        from which the structure is duplicated. In general the energy of the 'value' is lower than the
                        'key'

        """

        ret = {}
        selection = self.ids_sorted(ids)
        values = np.array([self.value(i) for i in selection])
        if len(values) == 0:
            return ret
        diffs = np.ediff1d(values)

        for i in range(len(diffs)):
            idiff = diffs[i]
            if idiff < self.value_tol:
                ident1 = selection[i]
                ident2 = selection[i + 1]
                pcm_log.debug('Testing distances between %s and %s' % (str(ident1), str(ident2)))
                distance = self.distance(ident1, ident2)
                # print 'Distance = ', distance
                if distance < self.distance_tolerance:
                    pcm_log.debug('Distance %7.3f < %7.3f' % (distance, self.distance_tolerance))
                    ret[ident2] = ident1
        if len(ret) > 0:
            pcm_log.debug('Number of duplicates %d' % len(ret))
        return ret

    def get_duplicates(self, ids, tolerance, fast=False):
        dupes_dict = {}
        dupes_list = []
        values = {}
        for i in ids:
            values[i] = self.value(i)
        selection = self.ids_sorted(ids)
        print('Searching duplicates in %d structures' % len(selection))
        for i in range(len(selection) - 1):
            entry_id = selection[i]
            value_i = values[entry_id]
            for j in range(i + 1, len(selection)):
                entry_jd = selection[j]
                if fast and entry_jd in dupes_list:
                    continue
                value_j = values[entry_jd]
                if abs(value_i - value_j) < self.value_tol:
                    distance = self.distance(entry_id, entry_jd)
                    if distance < tolerance:
                        if entry_id in dupes_dict:
                            dupes_dict[entry_id].append(entry_jd)
                        else:
                            dupes_dict[entry_id] = [entry_jd]
                        dupes_list.append(entry_jd)
        return dupes_dict, [x for x in selection if x in dupes_list]

    def cleaned_from_duplicates(self, ids):
        selection = self.ids_sorted(ids)
        duplicates_dict = self.check_duplicates(selection)
        return [x for x in selection if x not in duplicates_dict.keys()]

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

        ids_pair = [entry_id, entry_jd]
        ids_pair.sort()
        distance_entry = self.pcdb.db.distances.find_one({'pair': ids_pair}, {'distance': 1})
        self.pcdb.db.distances.create_index([("pair", ASCENDING)])

        if distance_entry is None:
            print('Distance not in DB')
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
            distance = float(np.mean(dij))
            self.pcdb.db.distances.insert({'pair': ids_pair, 'distance': distance})
        else:
            distance = distance_entry['distance']
        return distance

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
                print('Adding entry ' + str(entry['_id']) + ' from ' + readdb.name)
                self.new_entry(readdb.get_structure(entry['_id']))
                index += 1

    def move_random(self, entry_id, factor=0.2, in_place=False, kind='move'):
        structure = self.get_structure(entry_id)
        changer = StructureChanger(structure=structure)
        if kind == 'move':
            changer.random_move_many_atoms(epsilon=factor)
        else:  # change
            changer.random_change(factor)
        if in_place:
            return self.set_structure(entry_id, changer.new_structure)
        else:
            return self.new_entry(changer.new_structure, active=False)

    def move(self, entry_id, entry_jd, factor=0.2, in_place=False):
        """
        Moves entry_id in the direction of entry_jd
        If in_place is True the movement occurs on the
        same address as entry_id

        :param factor:
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
        match.match_size()
        match.match_shape()
        match.match_atoms()
        displacements = match.reduced_displacement()

        new_reduced = match.structure2.reduced + factor * displacements
        new_cell = match.structure2.cell
        new_symbols = match.structure2.symbols
        new_structure = Structure(reduced=new_reduced, symbols=new_symbols, cell=new_cell)
        if in_place:
            return self.set_structure(entry_id, new_structure)
        else:
            return self.new_entry(new_structure, active=False)

    def __str__(self):
        ret = '\n'
        ret += '[%s] Population type:    %s\n' % (self.tag, 'Relax Structures')
        ret += '[%s] Database:           %s\n' % (self.tag, self.name)
        ret += '[%s] Tag:                %s\n' % (self.tag, self.tag)
        ret += '[%s] Target-Forces:      %7.2E\n' % (self.tag, self.target_forces)
        ret += '[%s] Value tolerance:    %7.2E\n' % (self.tag, self.value_tol)
        ret += '[%s] Distance tolerance: %7.2E\n\n' % (self.tag, self.distance_tolerance)
        if self.composition is not None:
            ret += '[%s] Composition:                    %s\n' % (self.tag, self.composition.formula)
            ret += '[%s] Minimal composition multiplier: %d\n' % (self.tag, self.min_comp_mult)
            ret += '[%s] Maximal composition multiplier: %d\n' % (self.tag, self.max_comp_mult)
        ret += '[%s] Members:            %d\n' % (self.tag, len(self.members))
        ret += '[%s] Actives:            %d\n' % (self.tag, len(self.actives))
        ret += '[%s] Evaluated:          %d\n' % (self.tag, len(self.evaluated))
        ret += '\n'
        return ret

    def value(self, entry_id):
        entry = self.get_entry(entry_id)
        structure = self.get_structure(entry_id)
        if 'properties' not in entry:
            pcm_log.debug('This entry has no properties %s' % str(entry['_id']))
            return None
        elif entry['properties'] is None:
            return None
        elif 'energy' not in entry['properties']:
            pcm_log.debug('This entry has no energy in properties %s' % str(entry['_id']))
            return None
        else:
            return entry['properties']['energy'] / structure.get_composition().gcd

    @property
    def to_dict(self):
        return {'name': self.name,
                'tag': self.tag,
                'target_forces': self.target_forces,
                'value_tol': self.value_tol,
                'distance_tolerance': self.distance_tolerance}

    def from_dict(self, population_dict):
        return RelaxStructures(name=population_dict['name'],
                               tag=population_dict['tag'],
                               target_forces=population_dict['target_forces'],
                               value_tol=population_dict['value_tol'],
                               distance_tolerance=population_dict['distance_tolerance'])

    def cross(self, ids):

        assert len(ids) == 2

        structure1 = self.get_structure(ids[0])
        structure2 = self.get_structure(ids[1])

        split_match = SplitMatch(structure1, structure2)
        st1, st2 = split_match.get_simple_match()

        entry_id = self.new_entry(st1, active=True)
        entry_jd = self.new_entry(st2, active=True)

        return entry_id, entry_jd

    def str_entry(self, entry_id):

        struct = self.get_structure(entry_id)
        return str(struct)
