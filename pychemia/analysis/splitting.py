import itertools

import numpy as np

from pychemia import Composition, Structure


class SplitMatch:
    def __init__(self, structure1, structure2):

        assert structure1.nsites == structure2.nsites

        self.structures = [None, None]
        self.structures[0] = structure1
        self.structures[1] = structure2
        self.cut_planes = [None, None]
        self.split_sites = [None, None]
        self.ss_tags = [None, None]

    def get_minimal_splitting(self):

        delta_distance = 0.1
        while True:
            for i in range(2):
                self.cut_planes[i] = get_cut_planes(self.structures[i], delta_distance=delta_distance)
                self.split_sites[i] = get_simple_split_sites(self.structures[i], self.cut_planes[i])
                if len(self.split_sites[i]) == 0:
                    delta_distance *= 0.9
                    print('delta distance reduced to', delta_distance)
            if len(self.split_sites[0]) > 0 and len(self.split_sites[1]) > 0:
                # Select two compatible splits (that preserve the number of atoms)
                has_groups = False
                self.ss_tags[0] = None
                self.ss_tags[1] = None
                grp0 = None
                grp1 = None
                trial = 0
                while not has_groups:
                    rnd = np.random.randint(len(self.split_sites[0]))
                    grp0 = self.split_sites[0].keys()[rnd]
                    for grp1 in self.split_sites[1]:
                        trial += 1
                        print('Selected Groups:', grp0, grp1)
                        if len(self.split_sites[0][grp0]) == len(self.split_sites[1][grp1]):
                            has_groups = True
                            break
                    if trial > 10:
                        print('Bad splitting of sites, decreasing distance')
                        grp0 = None
                        grp1 = None
                        break

                if grp0 is not None and grp1 is not None:
                    self.ss_tags[0] = grp0
                    self.ss_tags[1] = grp1
                    break
                else:
                    delta_distance *= 0.9

    def get_simple_match(self):

        self.get_minimal_splitting()

        # grp1 and grp2 are tuples where the first value is the dimension of cut and the second one
        # is the second one is the index for the plane cut
        grp0 = self.ss_tags[0]
        grp1 = self.ss_tags[1]

        idim1 = grp0[0]
        idim2 = grp1[0]
        cp1 = grp0[1]
        cp2 = grp1[1]

        cut_value1 = self.cut_planes[0][idim1][cp1]
        cut_value2 = self.cut_planes[1][idim2][cp2]

        print('idim1:', idim1, 'cut_value1:', cut_value1)
        print('idim2:', idim2, 'cut_value2:', cut_value2)

        # Select one possible permutation that will make the second
        # structure being cutted along the same dimension as the first one
        permsel1 = None
        permsel2 = None
        for i in itertools.permutations(range(3), 3):
            if i[idim1] == idim2:
                permsel1 = i
            if i[idim2] == idim1:
                permsel2 = i
            if permsel1 is not None and permsel2 is not None:
                break

        print('Permutation Selected 1:', permsel1)
        print('Permutation Selected 2:', permsel2)

        newreduced1 = np.zeros((self.structures[0].natom, 3))
        newreduced2 = np.zeros((self.structures[1].natom, 3))
        symbols1 = self.structures[0].natom * ['U']
        symbols2 = self.structures[1].natom * ['U']

        nsites = self.structures[0].nsites

        index1 = 0
        index2 = 0

        for i in range(nsites):
            if i in self.split_sites[0][grp0]:
                newreduced1[index1] = self.structures[0].reduced[i, :]
                symbols1[index1] = self.structures[0].symbols[i]
                index1 += 1
            else:
                pos = self.structures[0].reduced[i, permsel1]
                # print i,' - ', self.structures[0].reduced[i], '==>', pos
                newreduced2[index2] = pos
                symbols2[index2] = self.structures[0].symbols[i]
                index2 += 1

        # print '1 Reduced\n', newreduced1
        # print '2 Reduced\n', newreduced2

        # print '1 Symbols:', symbols1
        # print '2 Symbols:', symbols2

        for i in range(nsites):
            if i in self.split_sites[1][grp1]:
                newreduced2[index2] = self.structures[1].reduced[i, :]
                symbols2[index2] = self.structures[1].symbols[i]
                index2 += 1
            else:
                pos = self.structures[1].reduced[i, permsel1]
                # print i,' - ', self.structures[1].reduced[i], '==>', pos
                newreduced1[index1] = pos
                symbols1[index1] = self.structures[1].symbols[i]
                index1 += 1

        cell1 = np.array(self.structures[0].cell)
        cell2 = np.array(self.structures[1].cell)

        # The Splitting is not checking for right stochiometry right now
        # Reverting the symbols to the original ones
        symbols1 = self.structures[0].symbols
        symbols2 = self.structures[1].symbols

        newst1 = Structure(symbols=symbols1, reduced=newreduced1, cell=cell1)
        newst2 = Structure(symbols=symbols2, reduced=newreduced2, cell=cell2)

        return newst1, newst2


def get_cut_planes(structure, delta_distance=0.1):
    """
    Compute all the possible cutting planes defined on reduced coordinates
    for which the structure can be efectively cut considering that the
    perpendicular distance between two atoms is larger than  'delta_distance'

    :param structure:
    :param delta_distance:
    :return:
    """
    ret = {}

    for idim in range(3):
        ret[idim] = []
        # Make a copy of the reduced positions
        v = np.array(structure.reduced[:, idim])
        # Extent the last position with a replica of the first one
        # Because the structure is periodic
        v.sort()
        v = np.concatenate((v, [1 + v[0]]))
        # Compute the vector of differences df[i]=v[i+1]-v[i]
        df = np.diff(v)
        for i in (range(structure.nsites)):
            if df[i] > delta_distance:
                ret[idim].append((v[i] + 0.5 * df[i]) % 1.0)
        ret[idim].sort()
    return ret


def get_simple_split_sites(structure, cut_planes):
    ret = {}
    for idim in range(3):
        for iplane in range(len(cut_planes[idim])):
            ret[(idim, iplane)] = []
            for i in range(structure.natom):
                if structure.reduced[i, idim] < cut_planes[idim][iplane]:
                    ret[(idim, iplane)].append(i)

    # Removing trivial partitions of the entire structure
    for i in ret.keys():
        if len(ret[i]) == structure.nsites:
            ret.pop(i)
    return ret


def get_split_sites(structure, idim, cut_planes, indices):
    """
    Return the indices of the sites corresponding to each
    partition of the structure

    :param structure:
    :param idim:
    :param cut_planes:
    :param indices:
    :return:
    """

    value0 = cut_planes[idim][indices[0]]
    value1 = cut_planes[idim][indices[1]]

    vrange = [[], []]
    if value1 > value0:
        vrange[0] = [(0, value0), (value1, 1)]
        vrange[1] = [(value0, value1)]
    else:
        vrange[0] = [(0, value1), (value0, 1)]
        vrange[1] = [(value1, value0)]

    factor0 = vrange[1][0][1] - vrange[1][0][0]
    factor1 = 1.0 - vrange[0][1][0] + vrange[0][0][1]

    ret = [[], []]
    for i in range(structure.nsites):
        for j in range(2):
            for irange in vrange[j]:
                if irange[0] <= structure.reduced[i, idim] < irange[1]:
                    ret[j].append(i)
    return ret, (factor0, factor1)


def get_all_splitted_compositions(structure, cut_planes):
    ret = []
    for idim in range(3):
        nplanes = len(cut_planes[idim])
        if nplanes > 1:
            for j in itertools.combinations(range(nplanes), 2):
                split_sites, factors = get_split_sites(structure, idim, cut_planes, j)
                npsymbols = np.array(structure.symbols)
                for k in range(2):
                    symbols = npsymbols[split_sites[k]]
                    comp = Composition(symbols)
                    factor = factors[k]
                    ret.append({'dim': idim, 'plane_indices': j, 'composition': comp.composition,
                                'natom': len(symbols), 'factor': factor})
    return ret


def get_matching_options(structure1, structure2, delta_distance=0.1):
    # Get Cutting Planes
    cp1 = get_cut_planes(structure1, delta_distance)
    cp2 = get_cut_planes(structure2, delta_distance)

    # Split Compositions
    sc1 = get_all_splitted_compositions(structure1, cp1)
    sc2 = get_all_splitted_compositions(structure2, cp2)

    # Indexing the Split Compositions
    index_sc1 = {}
    for i in range(len(sc1)):
        natom = sc1[i]['natom']
        if natom in index_sc1:
            index_sc1[natom].append(i)
        else:
            index_sc1[natom] = [i]

    index_sc2 = {}
    for i in range(len(sc2)):
        natom = sc2[i]['natom']
        if natom in index_sc2:
            index_sc2[natom].append(i)
        else:
            index_sc2[natom] = [i]

    ret = []
    for i in range(structure1.nsites / 2 + 1):
        if i in index_sc1 and i in index_sc2:
            for j in index_sc1[i]:
                for k in index_sc2[i]:
                    if sc1[j]['composition'] == sc2[k]['composition']:
                        ret.append({'dim1': sc1[j]['dim'], 'plane_indices1': sc1[j]['plane_indices'],
                                    'dim2': sc2[k]['dim'], 'plane_indices2': sc2[k]['plane_indices']})
    return ret, cp1, cp2


def do_mixing(structure1, structure2, cut_planes1, cut_planes2, matching, match_index):
    # Getting the atoms involved in both partitions:
    dim1 = matching[match_index]['dim1']
    partition1 = get_split_sites(structure1, dim1, cut_planes1, matching[match_index]['plane_indices1'])
    dim2 = matching[match_index]['dim2']
    partition2 = get_split_sites(structure1, dim2, cut_planes2, matching[match_index]['plane_indices2'])

    index1 = matching[match_index]['plane_indices1'][0]
    index2 = matching[match_index]['plane_indices1'][1]

    reduced11 = np.array([structure1.reduced[x] for x in partition1[0]])
    reduced11[:, dim1] = (reduced11[:, dim1] - cut_planes1[dim1][index1] + 1.0) % 1.0
    reduced11[:, dim1] -= np.min(reduced11[:, dim1])
    symbols = [structure1.symbols[x] for x in partition1[0]]
    st11 = Structure(reduced=reduced11, cell=structure1.cell, symbols=symbols)

    reduced12 = np.array([structure1.reduced[x] for x in partition1[1]])
    reduced12[:, dim1] = (reduced12[:, dim1] - cut_planes1[dim1][index2] + 1.0) % 1.0
    reduced12[:, dim1] -= np.min(reduced12[:, dim1])
    symbols = [structure1.symbols[x] for x in partition1[1]]
    st12 = Structure(reduced=reduced12, cell=structure1.cell, symbols=symbols)

    index1 = matching[match_index]['plane_indices2'][0]
    index2 = matching[match_index]['plane_indices2'][1]

    reduced21 = np.array([structure2.reduced[x] for x in partition2[0]])
    reduced21[:, dim2] = (reduced21[:, dim2] - cut_planes2[dim2][index1] + 1.0) % 1.0
    reduced21[:, dim2] -= np.min(reduced21[:, dim2])
    symbols = [structure2.symbols[x] for x in partition2[0]]
    st21 = Structure(reduced=reduced21, cell=structure2.cell, symbols=symbols)

    reduced22 = np.array([structure2.reduced[x] for x in partition2[1]])
    reduced22[:, dim2] = (reduced22[:, dim2] - cut_planes2[dim2][index2] + 1.0) % 1.0
    reduced22[:, dim2] -= np.min(reduced22[:, dim2])
    symbols = [structure2.symbols[x] for x in partition2[1]]
    st22 = Structure(reduced=reduced22, cell=structure2.cell, symbols=symbols)

    return st11, st12, st21, st22
