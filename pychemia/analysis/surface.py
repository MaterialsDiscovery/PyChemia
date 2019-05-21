
import numpy as np
import pychemia
import itertools
import scipy.spatial
from scipy.spatial import qhull
from pychemia.utils.periodic import covalent_radius


# return [x, y, d], that ax + by = d, d = gcd(a, b)
def ext_gcd(a, b):
    v1 = [1, 0, a]
    v2 = [0, 1, b]

    if a > b:
        a, b = b, a
    while v1[2] > 0:
        q = v2[2] / v1[2]
        for i in range(0, len(v1)):
            v2[i] -= q * v1[i]
        v1, v2 = v2, v1
    return v2


def print_vector(v):
    print('v = ', end='')
    for i in range(3):
        print(v[i], end='')
    print("\n")


def rotate_along_indices(structure, h, k, l, layers, tol=1.e-5):
    cell = structure.cell
    a1 = np.array(cell[0])
    a2 = np.array(cell[1])
    a3 = np.array(cell[2])
    rcell = structure.lattice.reciprocal().cell
    b1 = np.array(rcell[0])
    b2 = np.array(rcell[1])
    b3 = np.array(rcell[2])

    # Solve equation pk + ql = 1 for p and q using extended_Euclidean algorithm

    v = ext_gcd(k, l)
    p = v[0]
    q = v[1]
    # print('p = ', p)
    # print('q = ', q)

    k1 = np.dot(p * (k * a1 - h * a2) + q * (l * a1 - h * a3), l * a2 - k * a3)
    k2 = np.dot(l * (k * a1 - h * a2) - k * (l * a1 - h * a3), l * a2 - k * a3)
    # print("\n\nk1 = ", k1)
    # print("k2 = ", k2)

    if abs(k2) > tol:
        c = -int(round(k1 / k2))
        # print("c = -int(round(k1/k2)) = ", c)
        p, q = p + c * l, q - c * k

    # Calculate lattice vectors {v1, v2, v3} defining basis of the new cell

    v1 = p * (k * a1 - h * a2) + q * (l * a1 - h * a3)
    v2 = l * a2 - k * a3
    n = p * k + q * l
    v = ext_gcd(n, h)
    a = v[0]
    b = v[1]
    v3 = b * a1 + a * p * a2 + a * q * a3
    newbasis = np.array([v1, v2, v3])

    #    transformation = np.array([[p*k+q*l, -h*p, -h], [0, l, -k], [b, a*p, a*q]])
    #    inv_transformation = np.linalg.inv(transformation)

    symbols = []
    positions = structure.positions + np.tile(np.dot(structure.cell.T, np.array([0, 0, 0])).reshape(1, 3),
                                              (structure.natom, 1))
    reduced = np.linalg.solve(newbasis.T, positions.T).T
    for i in range(3):
        reduced[:, i] %= 1.0
    symbols += structure.symbols
    surf = pychemia.Structure(symbols=symbols, reduced=reduced, cell=newbasis)

    if layers > 1:
        new_surf = surf.supercell((1, 1, layers))
        cell = new_surf.cell
        # cell[2] = cell[2]+(layers+1)*cell[1]+(layers+1)*cell[0]
        surf = pychemia.Structure(symbols=new_surf.symbols, positions=new_surf.positions, cell=cell)

    #    print '\n\n********* Lattice vectors of the original cell *********\n\n', cell
    #    print '\n\n********* ATOMIC positions in the original cell **********\n', structure.positions
    #    print '\nTotal number of atoms in cell = ', structure.natom

    # Now create the surface starting from the original structure
    #    surf = structure.copy()
    #    surf.set_cell(newbasis)

    #    print '\n\n********* New basis of the surface cell *********\n\n', surf.cell
    #    print '\n\n********* Atomic coordinates in the newbasis of surface cell *********\n\n', surf.positions

    #    surf = surf.supercell((1, 1, layers))
    #    a1, a2, a3 = surf.cell

    #    surf.set_cell([a1, a2,
    #                   np.cross(a1, a2) * np.dot(a3, np.cross(a1, a2)) /
    #                   np.linalg.norm(np.cross(a1, a2)) ** 2])

    # Change unit cell to have the x-axis parallel with a surface vector
    # and z perpendicular to the surface:

    #    a1, a2, a3 = surf.cell
    #    surf.set_cell([(np.linalg.norm(a1), 0, 0),
    #                   (np.dot(a1, a2) / np.linalg.norm(a1),
    #                    np.sqrt(np.linalg.norm(a2) ** 2 - (np.dot(a1, a2) / np.linalg.norm(a1)) ** 2), 0),
    #                   (0, 0, np.linalg.norm(a3))])

    # Move atoms into the unit cell:

    #    scaled = surf.reduced
    #    scaled[:, :2] %= 1
    #    surf.set_reduced(scaled)
    #    surf.center(vacuum=vacuum, axis=2)
    return surf


def get_surface_atoms(structure):
    """
    Returns the list of atoms that belong to the surface of a given structure
    The surface atoms are computed as those atom for which the Voronoi tesselation is
    open, ie, their voronoi volume associated is infinite.

    :param structure: PyChemia Structure (Non-periodic in the current implementation)
    :return: (list) List of integers with the indices of the surface atoms
    """
    assert (not structure.is_periodic)
    voro = scipy.spatial.Voronoi(structure.positions)
    surface = [i for i in range(structure.natom) if -1 in voro.regions[voro.point_region[i]]]
    return surface


def get_surface_atoms_new(structure, use_covalent_radius=False):
    dln = scipy.spatial.Delaunay(structure.positions)

    if use_covalent_radius:
        simplices = []
        for j in dln.simplices:
            discard = False
            for ifacet in list(itertools.combinations(j, 3)):
                for ipair in itertools.combinations(ifacet, 2):
                    distance = np.linalg.norm(structure.positions[ipair[0]] - structure.positions[ipair[1]])
                    cov_distance = covalent_radius(structure.symbols[ipair[0]]) + covalent_radius(
                        structure.symbols[ipair[1]])
                    if distance > 3.0*cov_distance:
                        print('Distance: %f Cov-distance: %f' % (distance, cov_distance))
                        discard = True
                        break
            if not discard:
                print(j)
                simplices.append(j)
    else:
        simplices = dln.simplices

    c = np.array([[sorted(list(y)) for y in (itertools.combinations(x, 3))] for x in simplices])
    d = [list(x) for x in c.reshape((-1, 3))]

    ret = []
    dups = []
    for i in range(len(d) - 1):
        if d[i] in d[i+1:]:
            dups.append(d[i])
    for i in d:
        if i not in dups:
            ret.append(i)
    return np.unique(np.array(ret).flatten())


def get_onion_layers(structure):
    """
    Returns the different layers of a finite structure

    :param structure:
    :return:
    """
    assert (not structure.is_periodic)
    layers = []
    cur_st = structure.copy()

    morelayers = True
    while morelayers:
        pos = cur_st.positions
        if len(pos) <= 4:
            core = range(cur_st.natom)
            layers.append(core)
            break

        st = pychemia.Structure(positions=pos, symbols=len(pos)*['H'], periodicity=False)
        st.canonical_form()
        # print('The current volume is %7.3f' % st.volume)
        if st.volume < 0.1:
            core = range(cur_st.natom)
            layers.append(core)
            break

        try:
            voro = scipy.spatial.Voronoi(pos)
            surface = [i for i in range(cur_st.natom) if -1 in voro.regions[voro.point_region[i]]]
        except qhull.QhullError:
            surface = range(cur_st.natom)
            morelayers = False
        layers.append(surface)
        if not morelayers:
            break
        core = [i for i in range(cur_st.natom) if i not in surface]
        if len(core) == 0:
            break
        symbols = list(np.array(cur_st.symbols)[np.array(core)])
        positions = cur_st.positions[np.array(core)]
        cur_st = pychemia.Structure(symbols=symbols, positions=positions, periodicity=False)
    new_layers = [layers[0]]
    included = list(layers[0])
    acumulator = 0
    for i in range(1, len(layers)):
        noincluded = [j for j in range(structure.natom) if j not in included]
        # print 'Layer: %3d   Atoms on Surface: %3d     Internal Atoms: %3d' % (i,
        #                                                                      len(included) - acumulator,
        #                                                                      len(noincluded))
        acumulator = len(included)
        relabel = [noincluded[j] for j in layers[i]]
        included += relabel
        new_layers.append(relabel)

    return new_layers


def get_facets(structure, surface, seed, distance_tolerance=2.0):
    if seed not in surface:
        print('Error: seed not in surface')
        print('Seed: ', seed)
        print('Surface: ', surface)
        raise ValueError('Seed not in surface')

    idx_seed = surface.index(seed)
    pos = structure.positions[surface]
    dl = scipy.spatial.Delaunay(pos)

    # Delaunay tetrahedrals
    # The indices for each tetragon are relative to the surface
    # not the structure
    tetra = [x for x in dl.simplices if idx_seed in x]

    # Determining Facets
    facets = []
    for x in tetra:
        facets += [y for y in list(itertools.combinations(x, 3)) if idx_seed in y]

    # Indices of facets relative to surface not the structure
    facets = list(tuple([sorted(x) for x in facets]))

    selected_facets = []
    mintol = 1E5
    for x in facets:
        dm = scipy.spatial.distance_matrix(pos[x], pos[x])
        maxdist = max(dm.flatten())
        pair = np.where(dm == maxdist)
        atom1 = x[pair[0][0]]
        atom2 = x[pair[1][0]]
        covrad1 = pychemia.utils.periodic.covalent_radius(structure.symbols[surface[atom1]])
        covrad2 = pychemia.utils.periodic.covalent_radius(structure.symbols[surface[atom2]])
        if maxdist / (covrad1 + covrad2) < mintol:
            mintol = maxdist / (covrad1 + covrad2)
        if maxdist < distance_tolerance * (covrad1 + covrad2):
            # Converting indices from surface to structure
            st_facet = tuple([surface[x[i]] for i in range(3)])
            selected_facets.append(st_facet)
    return selected_facets, mintol + 0.1


def get_center_vector(structure, facet):
    v0 = np.array(structure.positions[facet[0]])
    v1 = np.array(structure.positions[facet[1]])
    v2 = np.array(structure.positions[facet[2]])
    facet_center = 1.0 / 3.0 * (v0 + v1 + v2)
    facet_vector = pychemia.utils.mathematics.unit_vector(np.cross(v1 - v0, v2 - v0))
    center_mass = structure.center_mass()
    if np.dot(facet_center - center_mass, facet_vector) < 0:
        facet_vector *= -1
    return (v0, v1, v2), facet_center, facet_vector


def attach_to_facet(structure, facet):
    (v0, v1, v2), facet_center, facet_vector = get_center_vector(structure, facet)

    rnd = np.random.random()
    if rnd < 0.4:
        return facet_center, facet_vector, tuple(facet)
    elif rnd < 0.6:
        edge_center = 0.5 * (v0 + v1)
        return edge_center, facet_vector, (facet[0], facet[1])
    elif rnd < 0.8:
        edge_center = 0.5 * (v1 + v2)
        return edge_center, facet_vector, (facet[1], facet[2])
    else:
        edge_center = 0.5 * (v0 + v2)
        return edge_center, facet_vector, (facet[0], facet[2])


def random_attaching(structure, seed, target_species, natom_crystal, radius=1.8, basetol=4.0):

    lys = pychemia.analysis.surface.get_onion_layers(structure)
    surface = lys[0]

    counter = 0
    while True:
        if seed not in surface or counter > 0:
            print('Current Seed not in surface, searching a new seed')
            seed = find_new_seed(structure, surface, seed, natom_crystal)

        tol = basetol
        facets = []
        while True:
            facets, mintol = get_facets(structure, surface, seed, distance_tolerance=tol)
            if len(facets) > 0:
                print('Possible Facets', facets)
                break
            elif mintol > 2 * basetol:
                return None, None, None, None
            else:
                tol = mintol
                print('No facets found, increasing tolerance to ', tol)

        counter = 0
        while True:
            counter += 1
            rnd = np.random.randint(len(facets))
            facet_chosen = facets[rnd]
            print('Seed: %3d     Number of facets: %3d     Facet chosen: %s' % (seed, len(facets), facet_chosen))

            center, uvector, atoms_facet = attach_to_facet(structure, facet_chosen)

            new_sts = {}
            good_pos = 0
            for specie in target_species:
                cov_rad = pychemia.utils.periodic.covalent_radius(specie)
                vec = center + radius * cov_rad * uvector
                new_symbols = list(structure.symbols) + [specie]
                new_positions = np.concatenate((structure.positions, [vec]))

                dist_matrix = scipy.spatial.distance_matrix(new_positions, new_positions)
                identity = np.eye(len(new_positions), len(new_positions))
                mindist = np.min((dist_matrix + 100 * identity).flatten())
                if mindist > cov_rad:
                    good_pos += 1
                    print('We have a minimal distance of', mindist)
                new_sts[specie] = pychemia.Structure(symbols=new_symbols, positions=new_positions, periodicity=False)

            if good_pos == len(target_species):
                print('Good position selected for all species')
                break
            else:
                print('No enough good positions: %d. One bad position, choosing a new facet' % good_pos)
                if counter > len(facets):
                    break
        if good_pos == len(target_species):
            break
    return new_sts, facet_chosen, center, uvector


def find_new_seed(st, surface, seed, natom_crystal):
    """
    Find a new atom to serve as seed for deposition

    :param st: Structure
    :param surface: Indices of atoms on surface
    :param seed: Current seed
    :param natom_crystal: Number of atoms in crystal, cluster atoms are always at the end
    :return:
    """

    candidates = [x for x in surface if x > natom_crystal]
    if seed in surface:
        new_seed = seed
    else:
        dists = scipy.spatial.distance_matrix(st.positions[seed].reshape((1, 3)),
                                              st.positions[surface])

        # First Option
        new_seed = surface[dists[0].argsort()[0]]

        # Second option
        dists = scipy.spatial.distance_matrix(st.center_mass().reshape((1, 3)),
                                              st.positions[candidates])

    if len(candidates) > 0:
        dists = scipy.spatial.distance_matrix(st.center_mass().reshape((1, 3)), st.positions[candidates])
        npcandidates = np.array(candidates)[dists[0].argsort()]
        print('Candidates ordered by distance to CM', npcandidates)
        for i in npcandidates:
            facets, mintol = get_facets(st, surface, i, distance_tolerance=1.0)
            if mintol < 5.0:
                new_seed = i
                break

    print('Seed: %3d => %3d' % (seed, new_seed))
    assert (new_seed in surface)
    return int(new_seed)
