__author__ = 'Guillermo Avendano-Franco'

from math import exp, sqrt
from sys import float_info
from itertools import combinations
from numpy import unique, array, zeros
from numpy.linalg import eig

from pychemia import Structure
from pychemia.db import StructureEntry
from pychemia.utils.periodic import atomic_number, covalent_radius, valence


class StructureAnalysis():
    """
    This class provides Structure Analysis.
    The set of analysis provided by this class uses only structural information of one
    single structure. The kind of analysis includes coordination numbers, bonds, distances,
    hardness, fingerprints.
    Most of those analysis rely on a robust computation of inter-atomic distances.
    This class uses lazy evaluation for the distances and bonds to lower cpu and memory footprint.
    """

    def __init__(self, structure):
        assert (isinstance(structure, Structure))
        self.structure = structure
        self._distances = None
        self._pairs = None

    def clean(self):
        """
        Delete the internal variables to allow them be recalculated again.

        :return: None
        """
        self._distances = None
        self._pairs = None

    def close_distances(self):
        """
        Computes the closest distances for all the atoms

        :return: (tuple) Return a bond's dictionary and distance's list
        """
        if self._pairs is None or self._distances is None:
            pairs_dict = {}
            distances_list = []
            index = 0
            for i, j in combinations(range(self.structure.natom), 2):
                ret = self.structure.get_cell().distance2(self.structure.reduced[i], self.structure.reduced[j])
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
                ret = self.structure.get_cell().distance2(self.structure.reduced[i], self.structure.reduced[i])
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
        
        #print self.structure.natom
        #print len(self._pairs)
        #print len(self._distances)
        return self._pairs, self._distances

    def get_bonds_coordination(self, tolerance=1.2, ensure_conectivity=False):
        """
        Computes simultaneously the bonds for all atoms and the coordination
        number using a multiplicative tolerance for the sum of covalent radius

        :param tolerance: (float) Tolerance factor (default is 1.2)
        :param ensure_conectivity: (bool) If True the tolerance of each bond is
               adjusted to ensure that each atom is connected at least once
        :return: tuple
        """
        bonds_dict, distances_list = self.close_distances()
        #print 'Bonds',len(bonds_dict)
        #print 'Dists',len(distances_list)
        bonds = []
        tolerances = []
        for i in range(self.structure.natom):
            tole = tolerance
            while True:
                tmp_bonds = []
                min_proportion = float_info.max
                for j in bonds_dict[str(i)]:
                    atom1 = self.structure.symbols[distances_list[j]['pair'][0]]
                    atom2 = self.structure.symbols[distances_list[j]['pair'][1]]
                    sum_covalent_radius = sum(covalent_radius([atom1, atom2]))
                    distance = distances_list[j]['distance']
                    if distance == 0.0:
                        continue
                    proportion = distance/sum_covalent_radius
                    min_proportion = min(min_proportion, proportion)
                    if proportion <= tole:
                        tmp_bonds.append(j)
                if len(tmp_bonds) == 0 and ensure_conectivity:
                    tole = min_proportion
                else:
                    bonds.append(tmp_bonds)
                    tolerances.append(min_proportion)
                    break

        coordination = [len(x) for x in bonds]
        return bonds, coordination, distances_list, tolerances

    def hardness(self, noupdate=False, verbose=False, tolerance=1.2):
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

        bonds, coordination, all_distances, tolerances = self.get_bonds_coordination(tolerance=tolerance,
                                                                                     ensure_conectivity=True)

        if verbose:
            print 'BONDS'
            print bonds
            print 'COORDINATION'
            print coordination

        sigma = 3.0
        c_hard = 1300.0
        x = 1.
        tot = 0.0
        f_d = 0.0
        f_n = 1.0
        atomicnumbers = atomic_number(self.structure.species)

        if verbose:
            print atomicnumbers

        for i in atomicnumbers:
            f_d += valence(i) / covalent_radius(i)
            f_n *= valence(i) / covalent_radius(i)
        if verbose:
            print 'fd', f_d
            print 'fn', f_n
            print atomicnumbers
        if f_d == 0:
            return 0.0
        f = 1.0 - (len(atomicnumbers) * f_n ** (1.0 / len(atomicnumbers)) / f_d) ** 2

        # Selection of different bonds
        diff_bonds = unique(array(reduce(lambda xx, y: xx+y, bonds)))
        for i in diff_bonds:
            i1 = all_distances[i]['pair'][0]
            i2 = all_distances[i]['pair'][1]

            ei = valence(self.structure.symbols[i1]) / covalent_radius(self.structure.symbols[i1])
            ej = valence(self.structure.symbols[i2]) / covalent_radius(self.structure.symbols[i2])
            #print 'bond ->', sqrt(ei * ej), (coordination[i1] * coordination[i2]), all_distances[i]['distance']

            sij = sqrt(ei * ej) / (coordination[i1] * coordination[i2]) / all_distances[i]['distance']
            num_i_j_bonds = len([j for j in diff_bonds if i1 in all_distances[j]['pair'] and
                                i2 in all_distances[j]['pair']])
            #print 'sij', sij
            #print 'num_i_j_bonds', num_i_j_bonds
            tot += num_i_j_bonds
            x *= sij
            #print 'x', x

        if verbose:
            print("V:", self.structure.volume)
            print("f:", f)
            print("x:", x)

        #print 'len_bonds', len(diff_bonds)
        hardness_value = c_hard / self.structure.volume * (len(diff_bonds) * x ** (1. / (len(diff_bonds)))) * exp(-sigma * f)

        if verbose:
            print hardness_value

        return round(hardness_value, 3)

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
            lap_m = zeros((n, n))
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
            eigen = eig(lap_m)[0]

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

    def hardness_OLD(self, noupdate=False, verbose=False, tolerance=0.05):
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
            print 'BONDS'
            print dis_dic
            print 'COORDINATION'
            print coord

        for i in dis_dic.keys():
            i1 = dis_dic[i][2][0]
            i2 = dis_dic[i][2][1]

            ei = valence(superc.symbols[i1]) / covalent_radius(superc.symbols[i1])
            ej = valence(superc.symbols[i2]) / covalent_radius(superc.symbols[i2])

            sij = sqrt(ei * ej) / (coord[i1] * coord[i2]) / dis_dic[i][0]

            tot += dis_dic[i][1]
            x *= sij * dis_dic[i][1]

        if verbose:
            print("V:", volume)
            print("f:", f)
            print("x:", x)

        hardness_value = c_hard / volume * (len(dis_dic) * x ** (1. / (len(dis_dic)))) * exp(-sigma * f)

        if verbose:
            print hardness_value

        return round(hardness_value, 3)


class EntryAnalysis():
    """
    This class provides a set of methods that uses not only the structural information but also
    tags, philogenetic info and/or computed properties.
    This class is created providing a StructureEntry instead of a simple Structure object
    This class uses lazy evaluation
    """

    def __init__(self, entry):
        assert (isinstance(entry, StructureEntry))
        self.entry = entry

