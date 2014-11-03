__author__ = 'Guillermo Avendano Franco'

import random
import itertools
import numpy as np
from pychemia.utils.mathematics import unit_vector


class StructureChanger():

    def __init__(self, structure):

        self.old_structure = structure
        self.new_structure = structure.copy()

        self.operations = []

    def permutator(self, pair):

        self.new_structure.symbols[pair[0]] = self.old_structure.symbols[pair[1]]
        self.new_structure.symbols[pair[1]] = self.old_structure.symbols[pair[0]]

        self.operations.append({'permutator': (pair[0], pair[1])})

    def random_permutator(self, forbidden_list=None):

        pairs = list(itertools.combinations(self.old_structure.species, 2))

        if forbidden_list is None:
            flist = []
        else:
            flist = forbidden_list

        for pair in pairs:
            for fpair in flist:
                if pair == tuple(fpair) or pair == tuple(reversed(fpair)):
                    pairs.remove(pair)

        pair = random.choice(pairs)

        indices0 = [i for i, x in enumerate(self.old_structure.symbols) if x == pair[0]]
        index0 = random.choice(indices0)

        indices1 = [i for i, x in enumerate(self.old_structure.symbols) if x == pair[1]]
        index1 = random.choice(indices1)

        self.permutator((index0, index1))

    def deform_cell(self, stress_eps):

        stress = np.eye(3) + np.diag(stress_eps[:3]) + \
                 np.array([[0.0, stress_eps[3], stress_eps[4]],
                           [stress_eps[3], 0.0, stress_eps[5]],
                           [stress_eps[4], stress_eps[5], 0.0]])

        self.new_structure.set_cell(np.dot(stress, self.old_structure.cell))

        self.operations.append({'deform_cell': stress_eps})

    def random_deform_cell(self, diag=True, nondiag=True, maxdelta=0.01):

        # 6 random numbers between -maxdelta and maxdelta
        stress_eps = np.random.rand(6)*2*maxdelta-maxdelta

        if diag:
            stress_eps[:3] = 0
        if nondiag:
            stress_eps[-3:] = 0

        self.deform_cell(stress_eps)

    def move_one_atom(self, index, vector):

        self.new_structure.positions[index] += vector

        self.operations.append({'move_one_atom': (index, vector)})

    def random_move_one_atom(self, epsilon):

        index = random.randint(0, len(self.old_structure))
        vector = unit_vector(np.random.rand(3)) * epsilon

        self.move_one_atom(index, vector)

    def random_move_many_atoms(self, epsilon=0.01, forbidden_indices=None, forbidden_species=None):

        if forbidden_indices is None:
            findices = []
        else:
            findices = forbidden_indices

        if forbidden_species is None:
            fspec = []
        else:
            fspec = forbidden_species

        for iatom in range(len(self.old_structure)):
            if iatom not in findices and self.new_structure.symbols[iatom] not in fspec:
                self.random_move_one_atom(epsilon)
