import itertools
import math
import sys

import numpy as np
import numpy.linalg
import functools

from pychemia import Structure, pcm_log
from pychemia.utils.mathematics import integral_gaussian
from pychemia.utils.periodic import atomic_number, covalent_radius, valence, atomic_symbol
from collections import OrderedDict


class StructureAnalysis:
    """
    The set of analysis provided by this class uses only structural information of one
    single structure. The kind of analysis includes coordination numbers, bonds, distances,
    hardness, fingerprints.
    Most of those analysis rely on a robust computation of inter-atomic distances.
    This class uses lazy evaluation for the distances and bonds to lower cpu and memory footprint.
    """

    def __init__(self, structure, supercell=(1, 1, 1), radius=50):
        """
        Takes one pychemia Structure object and will create a StructureAnalysis object. This object is computes the
        distances between all atoms in the unit cell with replicas extending up to a distance given by 'radius'.
        A supercell could be needed in cases where full connectivity of the entire crystal is needed.

        :param structure: A pychemia Structure object
        :param supercell: A supercell to be created from the Structure (Default=(1,1,1))
        :param radius: Maximal distance computed between two atoms
        """

        assert (isinstance(structure, Structure))
        if supercell != (1, 1, 1):
            self.structure = structure.supercell(supercell)
        else:
            self.structure = structure.copy()

        self._distances = None
        self._all_distances = None
        self._pairs = None
        self._supercell = supercell
        self._radius = radius
        # log.debug('Supercell : ' + str(self._supercell))
        # log.debug('Radius    : %7.2f' % self._radius)

    @property
    def radius(self):
        return self._radius

    @radius.setter
    def radius(self, value):
        assert (value > 0.0)
        if value != self._radius:
            self._distances = None
            self._pairs = None
            self._all_distances = None
        self._radius = value

    @property
    def distances(self):
        return self._distances

    def close_distances(self):
        """
        Computes the closest distances for all the atoms

        :return: (tuple) Return a bond's dictionary and distance's list
        """
        if self._pairs is None or self._distances is None:

            if self.structure.is_periodic:
                pcm_log.debug('Computing distances from scratch...')
                pairs_dict = {}
                distances_list = []
                index = 0
                for i, j in itertools.combinations(range(self.structure.natom), 2):
                    if index % 100 == 0:
                        pcm_log.debug('Computing distance between atoms %d and %d' % (i, j))
                    ret = self.structure.lattice.distance2(self.structure.reduced[i], self.structure.reduced[j],
                                                           radius=self.radius)
                    for k in ret:
                        if str(i) not in pairs_dict:
                            pairs_dict[str(i)] = [index]
                        else:
                            pairs_dict[str(i)].append(index)
                        if str(j) not in pairs_dict:
                            pairs_dict[str(j)] = [index]
                        else:
                            pairs_dict[str(j)].append(index)
                        ret[k]['pair'] = (i, j)
                        distances_list.append(ret[k])
                        index += 1
                for i in range(self.structure.natom):
                    ret = self.structure.lattice.distance2(self.structure.reduced[i], self.structure.reduced[i])
                    for k in ret:
                        if str(i) not in pairs_dict:
                            pairs_dict[str(i)] = [index]
                        else:
                            pairs_dict[str(i)].append(index)
                        ret[k]['pair'] = (i, i)
                        distances_list.append(ret[k])
                        index += 1
                self._pairs = pairs_dict
                self._distances = distances_list
            else:
                dm = self.structure.distance_matrix()
                dm += np.eye(len(dm)) * max(dm.flatten())
                pairs_dict = {}
                distances_list = []
                for i in range(self.structure.natom):
                    index = dm[:, i].argmin()
                    pairs_dict[str(i)] = [index]
                    distances_list.append(dm[index, i])
                self._pairs = pairs_dict
                self._distances = distances_list

        return self._pairs, self._distances

    def all_distances(self):

        if self._all_distances is None:
            ret = {}
            if not self.structure.is_periodic:
                dist_matrix = self.structure.distance_matrix()
            for i, j in itertools.combinations_with_replacement(range(self.structure.natom), 2):
                pair = (i, j)
                if self.structure.is_periodic:
                    ret[pair] = self.structure.lattice.distances_in_sphere(self.structure.reduced[i],
                                                                           self.structure.reduced[j],
                                                                           radius=self.radius)
                else:
                    ret[pair] = dist_matrix[i, j]

            self._all_distances = ret

        return self._all_distances

    def all_distances_by_species(self):

        all_distances = self.all_distances()
        ret = OrderedDict()

        atom_numbers = atomic_number(self.structure.species)
        a = list(itertools.combinations_with_replacement(atom_numbers, 2))
        keys = sorted([tuple(sorted(list(x))) for x in a])
        for key in keys:
            ret[key] = []

        for ipair in all_distances:
            key = tuple(sorted(atomic_number([self.structure.symbols[ipair[0]], self.structure.symbols[ipair[1]]])))
            if self.structure.is_periodic:
                ret[key] = np.concatenate((ret[key], all_distances[ipair]['distance']))
            else:
                ret[key].append(all_distances[ipair])

        # Sorting arrays
        for key in ret:
            ret[key].sort()
            ret[key] = np.array(ret[key])

        return ret

    def structure_distances(self, delta=0.01, sigma=0.01, integrated=True):
        dist_spec = self.all_distances_by_species()
        discrete_rdf = OrderedDict()
        nbins = int((self.radius + 5 * delta) / delta)
        discrete_rdf_x = np.arange(0, nbins * delta, delta)
        for spec_pair in dist_spec:
            discrete_rdf[spec_pair] = np.zeros(nbins)
            positive_distances = dist_spec[spec_pair][dist_spec[spec_pair] > 0]
            # log.debug('Pair %s' % str(spec_pair))
            for Rij in positive_distances:
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

    def fp_oganov(self, delta=0.01, sigma=0.01):
        struc_dist_x, struc_dist = self.structure_distances(delta=delta, sigma=sigma)
        fp_oganov = struc_dist.copy()
        vol = self.structure.volume
        for spec_pair in struc_dist:
            for i in range(len(struc_dist[spec_pair])):
                specie0 = atomic_symbol(spec_pair[0])
                specie1 = atomic_symbol(spec_pair[1])
                number_atoms0 = self.structure.composition[specie0]
                number_atoms1 = self.structure.composition[specie1]
                fp_oganov[spec_pair][i] *= vol / (delta * number_atoms0 * number_atoms1)
                fp_oganov[spec_pair][i] -= 1
        return struc_dist_x, fp_oganov

    def bonds_coordination(self, initial_cutoff_radius=0.8, use_laplacian=True, jump=0.01, tol=1E-15):

        cutoff_radius = initial_cutoff_radius
        ad = self.all_distances()
        bonds = {}
        while True:
            laplacian = np.zeros((self.structure.natom, self.structure.natom), dtype=np.int8)

            for pair in ad:
                atom1 = self.structure.symbols[pair[0]]
                atom2 = self.structure.symbols[pair[1]]
                sum_covalent_radius = sum(covalent_radius([atom1, atom2]))
                condition = np.bitwise_and(ad[pair]['distance'] < cutoff_radius * sum_covalent_radius,
                                           ad[pair]['distance'] > 0)
                bonds[pair] = ad[pair]['distance'][condition]
                if len(bonds[pair]) > 0:
                    laplacian[pair[0], pair[1]] = -1
                    laplacian[pair[1], pair[0]] = -1
            for i in range(self.structure.natom):
                laplacian[i, i] = 0
                laplacian[i, i] = -sum(laplacian[i])

            if use_laplacian:
                if np.max(np.abs(laplacian)) == 0:
                    cutoff_radius += jump

                ev = numpy.linalg.eigvalsh(laplacian)

                if sum(ev < tol) > 1:
                    cutoff_radius += jump
                else:
                    break
            else:
                break
        coordination = np.zeros(self.structure.natom, dtype=int)
        for pair in bonds:
            coordination[pair[0]] += len(bonds[pair])
            coordination[pair[1]] += len(bonds[pair])

        return bonds, coordination, round(cutoff_radius, 3)

    def get_bonds_coordination(self, initial_cutoff_radius=0.8, ensure_conectivity=False, use_laplacian=True,
                               verbose=False, tol=1E-15, jump=0.01, use_jump=True):
        """
        Computes simultaneously the bonds for all atoms and the coordination
        number using a multiplicative tolerance for the sum of covalent radius

        :param use_jump:
        :param jump:
        :param tol:
        :param verbose:
        :param use_laplacian:
        :param initial_cutoff_radius: (float) Tolerance factor (default is 1.2)
        :param ensure_conectivity: (bool) If True the tolerance of each bond is
               adjusted to ensure that each atom is connected at least once
        :return: tuple
        """
        if verbose:
            print('Computing all distances...')
        bonds_dict, distances_list = self.close_distances()
        if verbose:
            print('Number of distances computed: ', len(distances_list))

        cutoff_radius = initial_cutoff_radius
        bonds = None
        coordination = None
        tolerances = None

        while True:
            if verbose:
                print('Current cutoff radius : ', cutoff_radius)
            bonds = []
            tolerances = []
            for i in range(self.structure.natom):
                tole = cutoff_radius
                while True:
                    tmp_bonds = []
                    min_proportion = sys.float_info.max
                    for j in bonds_dict[str(i)]:
                        atom1 = self.structure.symbols[distances_list[j]['pair'][0]]
                        atom2 = self.structure.symbols[distances_list[j]['pair'][1]]
                        sum_covalent_radius = sum(covalent_radius([atom1, atom2]))
                        distance = distances_list[j]['distance']
                        if distance == 0.0:
                            continue
                        proportion = distance / sum_covalent_radius
                        min_proportion = min(min_proportion, proportion)
                        if proportion <= tole:
                            tmp_bonds.append(j)
                    if len(tmp_bonds) == 0 and ensure_conectivity:
                        tole = min_proportion
                        cutoff_radius = tole
                    else:
                        bonds.append(tmp_bonds)
                        tolerances.append(min_proportion)
                        break

            if use_laplacian:
                size = (self.structure.natom, self.structure.natom)
                laplacian = np.zeros(size, dtype=np.int8)
                for listibonds in bonds:
                    for ibond in listibonds:
                        data = distances_list[ibond]
                        i = data['pair'][0]
                        j = data['pair'][1]
                        laplacian[i, j] = -1
                        laplacian[j, i] = -1
                        # print '%d %d' % (i,j)
                for i in range(self.structure.natom):
                    laplacian[i, i] = 0
                    laplacian[i, i] = -sum(laplacian[i])
                if verbose:
                    print(laplacian)
                if np.max(np.abs(laplacian)) == 0:
                    cutoff_radius += jump
                    if verbose:
                        print('The laplacian is all zero')
                        print('Increasing cutoff radius by ', jump, 'A\n')
                    continue

                # if verbose:
                # print laplacian
                # evals, evecs = scipy.sparse.linalg.eigsh(laplacian)
                ev = numpy.linalg.eigvalsh(laplacian)
                if verbose:
                    print('Number of Eigenvalues close to zero :', sum(ev < tol))
                    print('Lowest Eigenvalues :', ev)
                    # print 'Lowest Eigenvalues :', evals

                if sum(ev < tol) > 1 and use_jump:
                    cutoff_radius += jump
                    if verbose:
                        print('Increasing cutoff radius by ', jump, 'A\n')
                else:
                    increase = False
                    for i in bonds:
                        if sum(i) == 0 and use_jump:
                            increase = True
                    if increase:
                        cutoff_radius += jump
                        if verbose:
                            print('Increasing cutoff radius by', jump, 'A\n')
                    else:
                        break
            else:
                break

        if bonds is not None:
            coordination = [len(x) for x in bonds]
        return bonds, coordination, distances_list, tolerances, cutoff_radius

    def hardness_XX(self, initial_cutoff_radius=0.8, use_laplacian=True):
        """
        Implementation of Hardness algorithm:
        
        First-principles structural design of superhard materials
        J. Chem. Phys. 138, 114101 (2013); https://doi.org/10.1063/1.4794424
        Xinxin Zhang, et al.
        """

        bonds, coordination, cutoff_radius = self.bonds_coordination(initial_cutoff_radius=initial_cutoff_radius,
                                                                     use_laplacian=use_laplacian)

        sigma = 3.0
        c_hard = 1300.0
        xprod = 1.
        tot = 0.0
        f_d = 0.0
        f_n = 1.0
        atomicnumbers = atomic_number(self.structure.species)

        pcm_log.debug('Atomic numbers in the structure : %s' % str(atomicnumbers))

        for i in atomicnumbers:
            f_d += valence(i) / covalent_radius(i)
            f_n *= valence(i) / covalent_radius(i)

        if f_d == 0:
            pcm_log.debug('Returning zero as hardness. f_d= %10.3f' % f_d)
            return 0.0
        f = 1.0 - (self.structure.nspecies * f_n ** (1.0 / self.structure.nspecies) / f_d) ** 2

        diff_bonds = [x for x in bonds if len(bonds[x]) > 0]

        for pair in diff_bonds:
            i1 = pair[0]
            i2 = pair[1]

            ei = valence(self.structure.symbols[i1]) / covalent_radius(self.structure.symbols[i1])
            ej = valence(self.structure.symbols[i2]) / covalent_radius(self.structure.symbols[i2])

            for dij in bonds[pair]:
                sij = math.sqrt(ei * ej) / (coordination[i1] * coordination[i2]) / dij
                xprod *= sij

            num_i_j_bonds = len(bonds[pair])
            pcm_log.debug('Number of bonds for pair %s = %d' % (str(pair), num_i_j_bonds))
            tot += num_i_j_bonds

        vol = self.structure.volume

        pcm_log.debug("Structure volume: %7.3f" % vol)
        pcm_log.debug("Total number of bonds: %d" % tot)
        pcm_log.debug("Bonds: %s" % str(bonds))

        hardness_value = (c_hard / vol) * tot * (xprod ** (1. / tot)) * math.exp(-sigma * f)

        return round(hardness_value, 3), cutoff_radius, coordination

    def hardness(self, verbose=True, initial_cutoff_radius=0.8, ensure_conectivity=False, use_laplacian=True,
                 use_jump=True, tol=1E-15):
        """
        Calculates the hardness of a structure based in the model of XX
        We use the covalent radii from pychemia.utils.periodic.
        If noupdate=False
        the Laplacian matrix method is not used and rcut is 2*max(cov_radii)

        :param use_jump:
        :param ensure_conectivity:
        :param verbose: (bool) To print some debug info
        :param initial_cutoff_radius: (float)
        :param use_laplacian: (bool) If True, the Laplacian method is used
        :param tol: (float) Tolerance for considering two atoms bonded

        :rtype : (float)
        """
        if self._supercell == (1, 1, 1) and verbose:
            print('''Only internal connectivity can be ensure, for complete connectivity in the crystal you must use a
                  supercell at of (2,2,2)''')

        bonds, coordination, all_distances, tolerances, cutoff_radius = \
            self.get_bonds_coordination(initial_cutoff_radius=initial_cutoff_radius,
                                        ensure_conectivity=ensure_conectivity,
                                        use_laplacian=use_laplacian, verbose=verbose, use_jump=use_jump, tol=tol)

        if verbose:
            print('Structure coordination : ', coordination)

        sigma = 3.0
        c_hard = 1300.0
        x = 1.
        tot = 0.0
        f_d = 0.0
        f_n = 1.0
        atomicnumbers = atomic_number(self.structure.species)

        if verbose:
            print('Atomic numbers in the structure :', atomicnumbers)

        for i in atomicnumbers:
            f_d += valence(i) / covalent_radius(i)
            f_n *= valence(i) / covalent_radius(i)

        # if verbose:
        # print 'fd', f_d
        #    print 'fn', f_n
        #    print atomicnumbers
        if f_d == 0:
            return 0.0
        f = 1.0 - (len(atomicnumbers) * f_n ** (1.0 / len(atomicnumbers)) / f_d) ** 2

        # Selection of different bonds
        diff_bonds = np.unique(np.array(functools.reduce(lambda xx, y: xx + y, bonds)))
        if verbose:
            print('Number of different bonds : ', len(diff_bonds))

        for i in diff_bonds:
            i1 = all_distances[i]['pair'][0]
            i2 = all_distances[i]['pair'][1]

            ei = valence(self.structure.symbols[i1]) / covalent_radius(self.structure.symbols[i1])
            ej = valence(self.structure.symbols[i2]) / covalent_radius(self.structure.symbols[i2])
            # print 'bond ->', sqrt(ei * ej), (coordination[i1] * coordination[i2]), all_distances[i]['distance']

            sij = math.sqrt(ei * ej) / (coordination[i1] * coordination[i2]) / all_distances[i]['distance']
            num_i_j_bonds = len([j for j in diff_bonds if i1 in all_distances[j]['pair'] and
                                 i2 in all_distances[j]['pair']])
            # print 'sij', sij
            # print 'num_i_j_bonds', num_i_j_bonds
            tot += num_i_j_bonds
            x *= sij
            # print 'x', x

        vol = self.structure.volume
        if verbose:
            print("Structure volume:", vol)
            # print("f:", f)
            # print("x:", x)

        # print 'len_bonds', len(diff_bonds
        hardness_value = c_hard / vol * (len(diff_bonds) * x ** (1. / (len(diff_bonds)))) * math.exp(-sigma * f)

        return round(hardness_value, 3), cutoff_radius, coordination

    def get_bonds(self, radius, noupdate=False, verbose=False, tolerance=0.05):
        """
        Calculates bond lengths between all atoms within "radius".

        :param radius: (float) The radius of the sphere were the bonds will be computed
        :param noupdate: (bool) If the original radius should be increased until a connected structure is obtained
        :param verbose: (bool) Print some info about the number of bonds computed
        :param tolerance: (float) Tolerance for creation of a bond

        :rtype : dict The values of the bonds are stored in self.info['bonds']

        """
        # The number of atoms in the cell
        n = self.structure.natom
        coordination = None
        dis_dic = None

        if verbose:
            print('Testing with %s of atoms in cell' % n)
            print('Starting with R_cut = %s' % radius)

        while True:
            size = (n, n)
            lap_m = np.zeros(size)
            dis_dic = {}
            ndifbonds = 0
            found = False
            # lap_mat[i,j] = lap_mat[j,i]
            for i in range(n - 1):
                for j in range(i + 1, n):
                    dis = self.structure.get_distance(i, j, with_periodicity=True)
                    if dis < radius:
                        if len(dis_dic) != 0:
                            for kstr, kj in dis_dic.items():
                                tstr = "%s%s" % (self.structure.symbols[i],
                                                 self.structure.symbols[j])
                                if abs(kj[0] - dis) < tolerance and tstr in kstr:
                                    dis_dic[kstr][1] += 1
                                    found = True
                                    break
                        if not found:
                            ndifbonds += 1
                            kstr = "%s%s%s%s" % (self.structure.symbols[i],
                                                 self.structure.symbols[j],
                                                 self.structure.symbols[i],
                                                 ndifbonds)
                            dis_dic[kstr] = [dis, 1, [i, j]]

                        lap_m[i][j] = -1
                        lap_m[j][i] = -1

            # lap_map[i,i]
            coordination = []
            for i in range(n):
                lap_m[i, i] = abs(sum(lap_m[i]))
                coordination.append(lap_m[i, i])

            # Get eigenvalues and check how many zeros;
            eigen = numpy.linalg.eig(lap_m)[0]

            nzeros = 0
            for i in eigen:
                if round(i.real, 5) == 0.0:
                    nzeros += 1
            if nzeros == 1 or (noupdate and len(dis_dic) > 0):
                break
            else:
                radius += 0.02

            if radius > 10.0:
                print("Cut off radius > 10.0")
                break

        if verbose:
            print('New cut off = %s' % radius)
            print('Bonds:')
            for i, j in dis_dic.items():
                print("  %s %s" % (i, j))

        return radius, coordination, dis_dic

    def hardness_old(self, noupdate=False, verbose=False, tolerance=0.05):
        """
        Calculates the hardness of a structure based in the model of XX
        We use the covalent radii from pychemia.utils.periodic.
        If noupdate=False
        the Laplacian matrix method is not used and rcut is 2*max(cov_radii)

        :param noupdate: (bool) If True, the Laplacian method is used
        :param verbose: (bool) To print some debug info
        :param tolerance: (float)

        :rtype : (float)
        """

        superc = self.structure.copy()
        superc.supercell(2, 2, 2)
        structure_analisys = StructureAnalysis(superc)

        natom = superc.natom
        volume = superc.volume

        max_covalent_radius = max(covalent_radius(superc.symbols))
        if verbose:
            print('Number of atoms', natom)
            print('Volume         ', volume)
            print('Covalent rad max', max_covalent_radius)
        rcut, coord, dis_dic = structure_analisys.get_bonds(2.0 * max_covalent_radius, noupdate, verbose, tolerance)

        sigma = 3.0
        c_hard = 1300.0
        x = 1.
        tot = 0.0
        f_d = 0.0
        f_n = 1.0
        dic_atms = {}
        for i in superc.symbols:
            dic_atms[i] = atomic_number(i)

        for i in dic_atms.keys():
            f_d += valence(i) / covalent_radius(i)
            f_n *= valence(i) / covalent_radius(i)
        f = 1.0 - (len(dic_atms) * f_n ** (1.0 / len(dic_atms)) / f_d) ** 2

        if verbose:
            print('BONDS')
            print(dis_dic)
            print('COORDINATION')
            print(coord)

        for i in dis_dic.keys():
            i1 = dis_dic[i][2][0]
            i2 = dis_dic[i][2][1]

            ei = valence(superc.symbols[i1]) / covalent_radius(superc.symbols[i1])
            ej = valence(superc.symbols[i2]) / covalent_radius(superc.symbols[i2])

            sij = math.sqrt(ei * ej) / (coord[i1] * coord[i2]) / dis_dic[i][0]

            tot += dis_dic[i][1]
            x *= sij * dis_dic[i][1]

        if verbose:
            print("V:", volume)
            print("f:", f)
            print("x:", x)

        hardness_value = c_hard / volume * (len(dis_dic) * x ** (1. / (len(dis_dic)))) * math.exp(-sigma * f)

        if verbose:
            print(hardness_value)

        return round(hardness_value, 3)
