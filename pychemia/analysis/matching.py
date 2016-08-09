import numpy as np
from pychemia.utils.mathematics import lcm, shortest_triple_set
from pychemia import Structure
import itertools


class StructureMatch:
    def __init__(self, structure1, structure2):
        """
        Creates a structure match between 2 structures
        The structures will be change to match their number of
        atoms and the order of atoms inside such that the distances
        between atoms in equivalent positions on both structures
        is minimized.

        :param structure1: (Structure)
        :param structure2: (Structure)
        """
        assert (isinstance(structure1, Structure))
        assert (isinstance(structure2, Structure))
        assert structure1.is_perfect
        assert structure2.is_perfect
        self.structure1 = structure1.copy()
        self.structure2 = structure2.copy()
        self.base_lattice = self.structure1.lattice

    def match_size(self):

        assert self.structure1.is_crystal
        assert self.structure2.is_crystal

        gcd1 = self.structure1.get_composition().gcd
        gcd2 = self.structure2.get_composition().gcd

        sts = np.array(shortest_triple_set(lcm(gcd1, gcd2) / gcd1)).astype(int)
        supercell_multiples = sts[self.structure1.lattice.lengths.argsort()[::-1]]
        self.structure1 = self.structure1.supercell(supercell_multiples)

        sts = np.array(shortest_triple_set(lcm(gcd1, gcd2) / gcd2))
        supercell_multiples = sts[self.structure2.lattice.lengths.argsort()[::-1]]
        self.structure2 = self.structure2.supercell(supercell_multiples)

    def match_shape(self):

        self.structure1.canonical_form()
        self.structure2.canonical_form()
        assert (self.structure1.symbols == self.structure2.symbols)

    def match_atoms(self):

        if self.structure1.natom != self.structure2.natom:
            raise ValueError('Match the size first')

        best = {}
        for specie in self.structure1.species:
            selection = np.array(self.structure1.symbols) == specie

            distance_matrix, close_images = self.base_lattice.minimal_distances(self.structure1.reduced[selection],
                                                                                self.structure2.reduced[selection])

            min_trace = 1E10
            best[specie] = None
            if self.structure1.natom < 7:
                for i in itertools.permutations(range(len(distance_matrix))):
                    if distance_matrix[:, np.array(i)].trace() < min_trace:
                        min_trace = distance_matrix[:, np.array(i)].trace()
                        best[specie] = i
            else:
                # Only consider permutations of 2 positions
                if len(distance_matrix) > 1:
                    for ipar in itertools.permutations(range(len(distance_matrix)), 2):
                        i = list(range(len(distance_matrix)))
                        i[ipar[0]] = ipar[1]
                        i[ipar[1]] = ipar[0]
                        if distance_matrix[:, np.array(i)].trace() < min_trace:
                            min_trace = distance_matrix[:, np.array(i)].trace()
                            best[specie] = i
                    for ipar in itertools.permutations(range(len(distance_matrix)), 4):
                        i = list(range(len(distance_matrix)))
                        i[ipar[0]] = ipar[1]
                        i[ipar[1]] = ipar[0]
                        i[ipar[2]] = ipar[3]
                        i[ipar[3]] = ipar[2]
                        if distance_matrix[:, np.array(i)].trace() < min_trace:
                            min_trace = distance_matrix[:, np.array(i)].trace()
                            best[specie] = i
                else:
                    best[specie] = [0]

            print('For specie %s best permutation is %s' % (specie, str(best[specie])))

        best_permutation = np.zeros(self.structure1.natom, dtype=int)
        index = 0
        while index < self.structure1.natom:
            specie = self.structure1.symbols[index]
            selection = np.array(self.structure1.symbols) == specie
            best_permutation[selection] = index + np.array(best[specie])
            index += len(best[specie])

        self.structure2.sort_sites_using_list(best_permutation)

    def reduced_displacement(self):

        assert (self.structure1.symbols == self.structure2.symbols)
        assert (self.structure1.nsites == self.structure2.nsites)
        assert (self.structure1.natom == self.structure2.natom)
        ret = np.zeros((self.structure1.nsites, 3))

        distance_matrix, close_images = self.base_lattice.minimal_distances(self.structure1.reduced,
                                                                            self.structure2.reduced)

        for i in range(self.structure1.nsites):
            x1 = self.structure1.reduced[i]
            x2 = self.structure2.reduced[i] + close_images[i, i]
            ret[i] = x2 - x1

        return ret

    def cell_displacement(self):

        return np.dot(self.structure1.cell, np.linalg.inv(self.structure2.cell))

    def cartesian_distances(self):

        rd = self.reduced_displacement()
        ret = np.zeros(self.structure1.nsites)
        for i in range(self.structure1.nsites):
            ret[i] = np.dot(np.dot(rd[i], self.base_lattice.metric), rd[i])
            ret[i] = np.sqrt(ret[i])
        return ret
