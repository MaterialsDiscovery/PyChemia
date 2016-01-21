import numpy as np
import pychemia
import itertools
import scipy.spatial
from scipy.spatial import qhull


# return [x, y, d], that ax + by = d, d = gcd(a, b)
def ext_gcd(a, b):
    v1 = [1, 0, a]
    v2 = [0, 1, b]

    if a > b:
        a, b = b, a
    while v1[2] > 0:
        q = v2[2] / v1[2]
        for i in range(0, len(v1)):
            v2[i] = v2[i] - q * v1[i]
        v1, v2 = v2, v1
    return v2


def print_vector(v):
    print 'v = ',
    for i in range(3):
        print v[i],
    print "\n"


def create_surface(structure, h, k, l, layers, tol=1.e-5):
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
    print 'p = ', p
    print 'q = ', q

    k1 = np.dot(p * (k * a1 - h * a2) + q * (l * a1 - h * a3), l * a2 - k * a3)
    k2 = np.dot(l * (k * a1 - h * a2) - k * (l * a1 - h * a3), l * a2 - k * a3)
    print "\n\nk1 = ", k1
    print "k2 = ", k2

    if abs(k2) > tol:
        c = -int(round(k1 / k2))
        print "c = -int(round(k1/k2)) = ", c
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
            print 'No atoms in core'
            break
        symbols = list(np.array(cur_st.symbols)[np.array(core)])
        positions = cur_st.positions[np.array(core)]
        cur_st = pychemia.Structure(symbols=symbols, positions=positions, periodicity=False)

    new_layers = [layers[0]]
    included = list(layers[0])

    acumulator = 0
    for i in range(1, len(layers)):
        noincluded = [j for j in range(structure.natom) if j not in included]
        print 'Layer: %3d   Atoms on Surface: %3d     Internal Atoms: %3d' % (i,
                                                                              len(included) - acumulator,
                                                                              len(noincluded))
        acumulator = len(included)
        relabel = [noincluded[j] for j in layers[i]]
        included += relabel
        new_layers.append(relabel)

    return new_layers


def get_edges_and_facets(structure, surface, seed, distance_tolerance=2.0):
    if seed not in surface:
        print 'Error: seed not in surface'
        print 'Seed: ', seed
        print 'Surface: ', surface
        raise ValueError('Seed not in surface')

    idx_seed = surface.index(seed)
    pos = structure.positions[surface]
    dl = scipy.spatial.Delaunay(pos)
    if seed not in surface:
        print 'Error: Seed not on the surface'

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
    for x in facets:
        dm = scipy.spatial.distance_matrix(pos[x], pos[x])
        maxdist = max(dm.flatten())
        pair = np.where(dm == maxdist)
        atom1 = x[pair[0][0]]
        atom2 = x[pair[1][0]]
        covrad1 = pychemia.utils.periodic.covalent_radius(structure.symbols[surface[atom1]])
        covrad2 = pychemia.utils.periodic.covalent_radius(structure.symbols[surface[atom2]])
        if maxdist < distance_tolerance * (covrad1 + covrad2):
            # Converting indices from surface to structure
            st_facet = tuple([surface[x[i]] for i in range(3)])
            selected_facets.append(st_facet)
        else:
            pass
            # print 'Excluded: ', atom1, atom2, surface[atom1], surface[atom2]

    print 'Original facets', facets
    facets = selected_facets
    print 'Selected facets', facets

    # Deternmining Edges
    edges = []
    for x in facets:
        # Make indices are relative to structure
        edges += [y for y in list(itertools.permutations(x, 2)) if y[0] == seed]
    edges = list(set(edges))
    return edges, facets


def attach_to_edge(structure, edges, edge_index=None, radius=1.0):
    if edge_index is None:
        rnd = np.random.randint(len(edges))
        chosen = edges[rnd]
    else:
        assert edge_index < len(edges)
        chosen = edges[edge_index]

    print 'Edge chosen: ', chosen
    v1 = structure.positions[chosen[0]]
    v2 = structure.positions[chosen[1]]
    vector = 0.5 * (v2 + v1)
    dv = 0.5 * pychemia.utils.mathematics.unit_vector(vector)
    n = 0
    cur_vector = None
    increase = True
    while True:
        cur_vector = vector + n * dv
        distances = scipy.spatial.distance_matrix(structure.positions, cur_vector.reshape((1, 3)))
        for j in range(structure.natom):
            sum_covalent_radius = pychemia.utils.periodic.covalent_radius(structure.symbols[j]) + radius
            if distances[j] < sum_covalent_radius:
                print 'Sum Cov Rad: %7.3f  Distance: %7.3f ' % (sum_covalent_radius, distances[j])
                n += 1
                increase = True
                break
            increase = False
        if not increase:
            break
    return cur_vector, chosen


def attach_to_facet(structure, facets, facet_index=None, radius=1.0):
    if facet_index is None:
        rnd = np.random.randint(len(facets))
        chosen = facets[rnd]
    else:
        assert facet_index < len(facets)
        chosen = facets[facet_index]

    print 'Facet chosen: ', chosen
    v1 = np.array(structure.positions[chosen[0]])
    v2 = np.array(structure.positions[chosen[1]])
    v3 = np.array(structure.positions[chosen[2]])
    center = 0.5 * (v3 + v2 + v1)
    dv = 0.5 * pychemia.utils.mathematics.unit_vector(np.cross(v2 - v1, v3 - v1))
    if np.dot(dv, center) < 0:
        dv = -dv
    n = 0
    increase = True
    cur_vector = None
    while True:
        cur_vector = center + n * dv
        distances = scipy.spatial.distance_matrix(structure.positions, cur_vector.reshape((1, 3)))
        for j in range(structure.natom):
            sum_covalent_radius = pychemia.utils.periodic.covalent_radius(structure.symbols[j]) + radius
            if distances[j] < sum_covalent_radius:
                print 'Sum Cov Rad: %7.3f  Distance: %7.3f ' % (sum_covalent_radius, distances[j])
                n += 1
                increase = True
                break
            increase = False
        if not increase:
            break
    return cur_vector, chosen


def random_attaching(structure, seed, target_species, radius=1.0):
    print 'Seed: ', seed

    lys = pychemia.analysis.surface.get_onion_layers(structure)
    surface = lys[0]

    tol = 2.0
    edges = []
    facets = []
    while True:
        edges, facets = get_edges_and_facets(structure, surface, seed, distance_tolerance=tol)
        if len(edges) > 0 and len(facets) > 0:
            break
        else:
            print 'No facets or edges found, increasing tolerance tol=', tol
            tol += 0.1

    print 'Number of Edges:', len(edges)
    print 'Number of Facets: ', len(facets)

    if True or np.random.randint(2) == 0:
        print 'Attaching to edge'
        vec, edge_facet = attach_to_edge(structure, edges, radius=radius)
        attach_to = 'Edge'
    else:
        print 'Attaching to facet'
        vec, edge_facet = attach_to_facet(structure, facets, radius=radius)
        attach_to = 'Facet'

    new_sts = {}
    for specie in target_species:
        new_symbols = list(structure.symbols) + [specie]
        new_positions = np.concatenate((structure.positions, [vec]))

        new_sts[specie] = pychemia.Structure(symbols=new_symbols, positions=new_positions, periodicity=False)
    return new_sts, attach_to, edge_facet
