import math
import numpy as np
import scipy.spatial
import itertools
from pychemia.utils.mathematics import integral_gaussian
from pychemia.utils.periodic import atomic_number

__author__ = "Guillermo Avendano-Franco"


class ClusterAnalysis:
    def __init__(self, structure):
        """
        Cluster Analysis provides routines to compute Structure Analysis
        for finite systems such as molecules and clusters.

        :param structure: (pychemia.Structure) A PyChemia Structure object

        """

        self.structure = structure.copy()
        self.structure.relocate_to_cm()
        # This is the maximal distance between any two atoms
        self.max_distance = np.max(self.distance_matrix())

    def distance_matrix(self):
        """
        Compute the distance matrix, the distance between any two particles
        in the entire structure.
        For distance matrices related separated by species, use instead
        distance_matrix_by_species.

        :return:
        """
        return scipy.spatial.distance_matrix(self.structure.positions, self.structure.positions)

    def distance_matrix_by_species(self, specie1, specie2):

        assert self.structure.is_perfect
        group1 = np.array(self.structure.symbols) == specie1
        group2 = np.array(self.structure.symbols) == specie2
        return scipy.spatial.distance_matrix(self.structure.positions[group1], self.structure.positions[group2])

    def all_distances_by_species(self):

        ret = {}
        for i, j in itertools.combinations_with_replacement(self.structure.species, 2):
            pair = tuple(np.sort([i, j]))
            dm = self.distance_matrix_by_species(i, j)
            distances = np.triu(dm, 1).flatten()
            distances = np.sort(distances[distances > 0.0])
            ret[pair] = distances
        return ret

    def discrete_radial_distribution_function(self, delta=0.01, sigma=0.01, integrated=True):
        dist_spec = self.all_distances_by_species()
        discrete_rdf = {}
        nbins = int((self.max_distance + 8 * sigma) / delta)
        if nbins > 10000:
            nbins = 10000
        discrete_rdf_x = np.arange(0, nbins * delta, delta)
        for spec_pair in dist_spec:
            discrete_rdf[spec_pair] = np.zeros(nbins)
            for Rij in dist_spec[spec_pair]:
                # Integrating for bin from - 8*sigma to +8*sigma centered on Rij
                # Values outside this range are negligible
                imin = int(max(0, (Rij - 8 * sigma) / delta))
                imax = int(min(len(discrete_rdf_x), (Rij + 8 * sigma) / delta))
                # log.debug('Computing for distance %7.3f for indices between %d and %d' % (Rij, imin, imax))
                for i in range(imin, imax):
                    x = discrete_rdf_x[i]
                    if not integrated:
                        discrete_rdf[spec_pair][i] += np.exp(-((x - Rij) ** 2) / (2 * sigma * sigma)) / (
                            4 * math.pi * Rij * Rij)
                    else:
                        discrete_rdf[spec_pair][i] += integral_gaussian(x, x + delta, Rij, sigma) / (
                            4 * math.pi * Rij * Rij)

        return discrete_rdf_x, discrete_rdf


class ClusterMatch:
    def __init__(self, structure1, structure2):
        """
        Given two structures, ClusterMatch can compute a one-to-one
        association between atoms of both structures, trying to minimize
        the displacements needed to move one structure into another.

        :param structure1:
        :param structure2:
        """

        assert not structure1.is_periodic
        assert not structure2.is_periodic

        self.structure1 = structure1.copy()
        self.structure2 = structure2.copy()

        self.structure1.canonical_form()
        self.structure2.canonical_form()

        assert self.structure1.symbols == self.structure2.symbols

    def match(self):

        species = self.structure1.species
        species = list(np.array(species)[np.array(atomic_number(species)).argsort()])

        permutation = np.zeros(len(self.structure1.positions), dtype=int)
        num = 0
        for ispecie in species:

            pos1 = self.structure1.positions[np.array(self.structure1.symbols) == ispecie]
            pos2 = self.structure2.positions[np.array(self.structure2.symbols) == ispecie]
            dm = scipy.spatial.distance_matrix(pos1, pos2)

            match_list = np.zeros(len(pos1))
            maxdis = np.max(dm.flatten())
            for i in range(len(pos1)):
                match = dm[i, :].argmin()
                # print " %3s %3d %9.3f %9.3f %9.3f" % (ispecie, match, np.linalg.norm(pos1[i]),
                # np.linalg.norm(pos1[match]), dm[i,match])
                match_list[i] = match
                dm[:, match] = maxdis + 1

            match_list += num
            permutation[num:num + len(pos1)] = match_list
            num += len(pos1)

        # print permutation
        self.structure2.sort_sites_using_list(permutation)
