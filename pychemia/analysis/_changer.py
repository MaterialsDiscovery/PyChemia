import random
import itertools
import numpy as np
from pychemia.utils.mathematics import apply_rotation


class StructureChanger:
    def __init__(self, structure):

        self.old_structure = structure
        self.new_structure = structure.copy()

        self.operations = []

    def permutator(self, pair):

        self.new_structure.symbols[pair[0]] = self.old_structure.symbols[pair[1]]
        self.new_structure.symbols[pair[1]] = self.old_structure.symbols[pair[0]]

        self.operations.append({'permutator': (pair[0], pair[1])})

    def random_permutator(self, forbidden_list=None):

        assert (self.old_structure.nspecies > 1)
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

        stress = np.eye(3) + np.diag(stress_eps[:3]) + np.array([[0.0, stress_eps[3], stress_eps[4]],
                                                                 [stress_eps[3], 0.0, stress_eps[5]],
                                                                 [stress_eps[4], stress_eps[5], 0.0]])

        self.new_structure.set_cell(np.dot(stress, self.old_structure.cell))

        self.operations.append({'deform_cell': stress_eps})

    def random_deform_cell(self, diag=True, nondiag=True, maxdelta=0.01):

        # 6 random numbers between -maxdelta and maxdelta
        stress_eps = np.random.random(6) * 2 * maxdelta - maxdelta

        if diag:
            stress_eps[:3] = 0
        if nondiag:
            stress_eps[-3:] = 0

        self.deform_cell(stress_eps)

    def move_one_atom(self, index, vector):

        self.new_structure.positions[index] += vector
        self.new_structure.positions2reduced()

        self.operations.append({'move_one_atom': (index, vector)})

    def random_move_one_atom(self, mu=0.1, sigma=0.01):

        index = random.randint(0, len(self.old_structure) - 1)
        radius = np.abs(np.random.normal(mu, sigma))
        theta_x = 2 * np.pi * np.random.random_sample()
        theta_y = 2 * np.pi * np.random.random_sample()
        theta_z = 2 * np.pi * np.random.random_sample()
        vector = apply_rotation([1, 0, 0], theta_x, theta_y, theta_z)

        self.move_one_atom(index, vector * radius)

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

    def random_change(self, epsilon):
        rnd = np.random.random()
        if rnd < 0.25:
            self.random_deform_cell(maxdelta=epsilon)
        elif rnd < 0.5:
            self.random_move_many_atoms(epsilon=epsilon)
        elif rnd < 0.75 and self.old_structure.nspecies > 1:
            self.random_permutator()
        else:
            self.random_move_one_atom(mu=epsilon, sigma=0.01)
