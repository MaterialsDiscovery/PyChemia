"""
Delaunay Reduction
"""

import numpy as _np


def get_reduced_bases(cell, tolerance=1e-5):
    """
    This is an implementation of Delaunay reduction.
    Some information is found in International table.

    :param cell: (numpy.ndarray) Cell parameters dimension 3x3
    :param tolerance: (float) Tolerance used to reduce basis
    :return: (numpy.ndarray) a new cell
    """
    return get_delaunay_reduction(cell, tolerance)


def get_delaunay_reduction(lattice, tolerance):
    """
    Return a reduced cell using the delanuay reduction

    :param lattice: (numpy.ndarray) the cell parameters
    :param tolerance: (float) Tolerance used to reduce basis
    :return: (numpy.ndarray) a new cell
    """
    extended_bases = _np.zeros((4, 3), dtype=float)
    extended_bases[:3, :] = lattice
    extended_bases[3] = -_np.sum(lattice, axis=0)

    i = 0
    for i in range(100):
        if reduce_bases(extended_bases, tolerance):
            break
    if i == 99:
        print("Delaunary reduction is failed.")

    shortest = get_shortest_bases_from_extented_bases(extended_bases, tolerance)

    return shortest


def reduce_bases(extended_bases, tolerance):
    """
    Reduces the basis with a given tolerance

    :param extended_bases:
    :param tolerance: (float) Tolerance used by reduction algorithm
    """
    metric = _np.dot(extended_bases, extended_bases.transpose())
    for i in range(4):
        for j in range(i + 1, 4):
            if metric[i][j] > tolerance:
                for k in range(4):
                    if (not k == i) and (not k == j):
                        extended_bases[k] += extended_bases[i]
                extended_bases[i] = -extended_bases[i]
                extended_bases[j] = extended_bases[j]
                return False

    # Reduction is completed.
    # All non diagonal elements of metric tensor is negative.
    return True


def get_shortest_bases_from_extented_bases(extended_bases, tolerance):
    # print 'get_shortest_bases_from_extented_bases',tolerance

    def mycmp(x, y):
        return cmp(_np.vdot(x, x), _np.vdot(y, y))

    basis = _np.zeros((7, 3), dtype=float)
    basis[:4] = extended_bases
    basis[4] = extended_bases[0] + extended_bases[1]
    basis[5] = extended_bases[1] + extended_bases[2]
    basis[6] = extended_bases[2] + extended_bases[0]
    # Sort bases by the lengthes (shorter is earlier)
    basis = sorted(basis, cmp=mycmp)

    # Choose shortest and linearly independent three bases
    # This algorithm may not be perfect.
    for i in range(7):
        for j in range(i + 1, 7):
            for k in range(j + 1, 7):
                if abs(_np.linalg.det([basis[i], basis[j], basis[k]])) > tolerance:
                    return _np.array([basis[i], basis[j], basis[k]])

    print("Delaunary reduction is failed.")
    return basis[:3]
