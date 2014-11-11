"""
Mathematical operations
"""

__author__ = 'Guillermo Avendano-Franco'

import math
import numpy as _np
import itertools as _it


def length_vector(v):
    """
    Returns the length of a vector 'v'
    Arbitrary number of dimensions

    :param v: list, numpy.ndarray
    :rtype : float

    Example:
    >>> length_vector([1,2,3])
    3.7416573867739413
    """
    return _np.linalg.norm(v)


def length_vectors(m):
    """
    Returns the lengths of several vectors
    arranged as rows in a MxN matrix

    :param m: numpy.ndarray
    :rtype : object

    Example:
    >>> length_vectors([[1,2,3], [4,5,6], [7,8,9], [1,0,0], [0,0,2]])
    array([  3.74165739,   8.77496439,  13.92838828,   1.        ,   2.        ])
    """
    return _np.apply_along_axis(_np.linalg.norm, 1, m)


def unit_vector(v):
    """
    Returns the unit vector of the vector.
    Arbitrary number of dimensions

    :param v: list, numpy.array
    :rtype : numpy.ndarray
    Example:
>>> a = unit_vector([1, 2, 3])
>>> a
    array([ 0.26726124,  0.53452248,  0.80178373])
>>> length_vector(a)
    1.0
    """
    return _np.array(v) / length_vector(_np.array(v, dtype=float))


def unit_vectors(m):
    """
    Returns the unit vectors of a set
    of vectors arranged as rows in MxN matrix

    :param m: numpy.ndarray
    :rtype : numpy.ndarray

    Example:
>>> b = unit_vectors([[1, 2, 3], [4, 5, 6], [7, 8, 9], [1, 0, 0], [0, 0, 2]])
>>> b
array([[ 0.26726124,  0.53452248,  0.80178373],
       [ 0.45584231,  0.56980288,  0.68376346],
       [ 0.50257071,  0.57436653,  0.64616234],
       [ 1.        ,  0.        ,  0.        ],
       [ 0.        ,  0.        ,  1.        ]])
>>> length_vectors(b)
array([ 1.,  1.,  1.,  1.,  1.])
    """
    return _np.divide(_np.array(m, dtype=float).T, length_vectors(_np.array(m, dtype=float))).T


def angle_vector(v1, v2, units='rad'):
    """
    Returns the angle in radians (default) or degrees
    between vectors 'v1' and 'v2'::

    :param v1: (list, numpy.ndarray)
    :param v2: (list, numpy.ndarray)
    :param units: (str) : 'rad' (default) Radians
                          'deg' Degrees
    :rtype : float

    Examples:
    >>> angle_vector([1, 0, 0], [0, 1, 0])
    1.5707963267948966
    >>> angle_vector([1, 0, 0], [1, 0, 0])
    0.0
    >>> angle_vector([1, 0, 0], [-1, 0, 0])
    3.1415926535897931
    >>> angle_vector([1, 0, 0], [0, 1, 0],units='deg')
    90.0
    >>> angle_vector([1, 0, 0], [-1, 0, 0], units='deg')
    180.0
    """
    assert(units in ['rad', 'deg'])

    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    angle = _np.arccos(_np.dot(v1_u, v2_u))
    if _np.isnan(angle):
        if (v1_u == v2_u).all():
            return 0.0
        else:
            return _np.pi
    if units == 'rad':
        return angle
    elif units == 'deg':
        return 180.0 * angle / _np.pi


def angle_vectors(m, units='rad'):
    """
    Returns all the angles for all the
    vectors arranged as rows in matrix 'm'

    :param m: (numpy.ndarray)
    :param units: (str) : 'rad' Radians
                          'deg' Degrees

    :rtype : numpy.ndarray
    Example:

>>> a = angle_vectors([[1, 2, 3], [4, 5, 6], [7, 8, 9], [1, 0, 0], [0, 0, 2]])
>>> import pprint
>>> pprint.pprint(a)
    {(0, 1): 0.22572612855273419,
     (0, 2): 0.2858867976945072,
     (0, 3): 1.3002465638163236,
     (0, 4): 0.6405223126794245,
     (1, 2): 0.060160669141772885,
     (1, 3): 1.0974779950809703,
     (1, 4): 0.8178885561654512,
     (2, 3): 1.0442265974045177,
     (2, 4): 0.86825103780276369,
     (3, 4): 1.5707963267948966}
>>> a = angle_vectors([[1, 2, 3], [4, 5, 6], [7, 8, 9], [1, 0, 0], [0, 0, 2]], units='deg')
>>> pprint.pprint(a)
    {(0, 1): 12.933154491899135,
     (0, 2): 16.380106926405656,
     (0, 3): 74.498640433063002,
     (0, 4): 36.699225200489877,
     (1, 2): 3.4469524345065143,
     (1, 3): 62.880857226618922,
     (1, 4): 46.861562380328941,
     (2, 3): 59.829776886585428,
     (2, 4): 49.747120023952057,
     (3, 4): 90.0}
    """

    ret = {}
    for i in _it.combinations(range(len(m)), 2):
        ret[i] = angle_vector(m[i[0]], m[i[1]], units=units)
    return ret


def distance(v1, v2):
    """
    Return the vector v2-v1, the vector going from v1 to v2
    and the magnitude of that vector.

    :param v1: (list, numpy.ndarray)
    :param v2: (list, numpy.ndarray)
    :rtype : tuple

    Example:
    >>> distance([0,0,0,1],[1,0,0,0])
    (array([ 1,  0,  0, -1]), 1.4142135623730951)
    >>> distance([-1,0,0],[1,0,0])
    (array([2, 0, 0]), 2.0)
    """
    ret = _np.array(v2)-_np.array(v1)
    return ret, length_vector(ret)


def distances(m):
    """
    Return all the distances for all possible combinations
    of the row vectors in matrix m

    :param m: (list, numpy.ndarray)
    :rtype : dict

    Example:

>>> import pprint
>>> pprint.pprint(distances([[1,2,3], [4,5,6], [7,8,9], [1,0,0], [0,0,2]]))
    {(0, 1): (array([3, 3, 3]), 5.196152422706632),
     (0, 2): (array([6, 6, 6]), 10.392304845413264),
     (0, 3): (array([ 0, -2, -3]), 3.6055512754639891),
     (0, 4): (array([-1, -2, -1]), 2.4494897427831779),
     (1, 2): (array([3, 3, 3]), 5.196152422706632),
     (1, 3): (array([-3, -5, -6]), 8.3666002653407556),
     (1, 4): (array([-4, -5, -4]), 7.5498344352707498),
     (2, 3): (array([-6, -8, -9]), 13.45362404707371),
     (2, 4): (array([-7, -8, -7]), 12.727922061357855),
     (3, 4): (array([-1,  0,  2]), 2.2360679774997898)}
    """
    ret = {}
    for i in _it.combinations(range(len(m)), 2):
        ret[i] = distance(m[i[0]], m[i[1]])
    return ret


def wrap2_pmhalf(x):
    """
    Wraps a number or array in the interval ]-1/2, 1/2]
    values = -1/2 will be wrapped  to 1/2

    :param x:

    Example:

    >>> wrap2_pmhalf(-0.5)
    0.5
    >>> wrap2_pmhalf(0.0)
    0.0
    >>> wrap2_pmhalf([-0.75,-0.5,-0.25,0.0,0.25,0.5,0.75])
    array([ 0.25,  0.5 , -0.25,  0.  ,  0.25,  0.5 , -0.25])
    >>> wrap2_pmhalf([[-0.75,-0.5,-0.25],[0.25,0.5,0.75]])
    array([[ 0.25,  0.5 , -0.25],
           [ 0.25,  0.5 , -0.25]])
    """

    def wrap(num):
        tol12 = 1e-12
        if num > 0:
            ret = (num+0.5-tol12) % 1.0 - 0.5 + tol12
        else:
            ret = -(-(num - 0.5 - tol12) % 1.0) + 0.5 + tol12
        for y in [-0.25, 0.0, 0.25, 0.5]:
            ret = (lambda num2: y if abs(y-num2) < tol12 else num2)(ret)
        return ret

    if _np.iterable(x):
        vec = _np.vectorize(wrap)
        return vec(x)
    else:
        return wrap(x)


def vector_set_perpendicular(vector3):
    """
    Produces a set of three mutually perpendicular vectors
    The two other vectors will be unitary

    :return: (tuple) Two numpy arrays
    """
    v1 = unit_vector(vector3)
    v2 = None
    v3 = None
    while True:
        other = unit_vector(_np.random.rand(3))
        if _np.abs(_np.dot(v1, other)) > 0.05:
            v2 = unit_vector(_np.cross(v1, other))
            v3 = unit_vector(_np.cross(v1, v2))
            break
        else:
            continue
    #print _np.dot(v1, v2)
    #print _np.dot(v1, v3)
    #print _np.dot(v2, v3)
    assert (_np.abs(_np.dot(v1, v2)) < 1E-15)
    assert (_np.abs(_np.dot(v1, v3)) < 1E-15)
    assert (_np.abs(_np.dot(v2, v3)) < 1E-15)
    return v1, v2, v3


def matrix_from_eig(v1, v2, v3, lam1, lam2, lam3):
    """
    Given 3 eigenvectors and 3 eigenvalues, returns the
    matrix A that has those eigenvectors and eigenvalues.

    The matrix $A = P.D.P^{-1}$

    Where P is the column stack of eigenvectors and
    D is a diagonal matrix of eigevalues

    :param v1: First eigenvector
    :param v2: Second eigenvector
    :param v3: Third eigenvector
    :param lam1: First eigenvalue
    :param lam2: Second eigenvalue
    :param lam3: Third eigenvalue
    :return: (numpy.ndarray) The matrix
    """
    matrixp = _np.vstack((v1, v2, v3)).T
    matrixd = _np.diag([lam1, lam2, lam3])
    matrixpinv = _np.linalg.inv(matrixp)
    matrixa = _np.dot(matrixp, _np.dot(matrixd, matrixpinv))
    return matrixa


def integral_gaussian(a, b, mu, sigma):
    """
    Computes the integral of a gaussian centered
    in mu with a given sigma
    :param a:
    :param b:
    :param mu:
    :param sigma:
    :return:
    """

    #Integral from -\infty to a
    val_floor = 0.5 * (1+math.erf((a-mu)/(sigma*math.sqrt(2.0))))

    #Integral from -\infty to b
    val_ceil = 0.5 * (1+math.erf((b-mu)/(sigma*math.sqrt(2.0))))

    return val_ceil - val_floor

def frexp10(x):
    exp = int(math.floor(math.log10(abs(x))))
    return x / 10**exp, exp

def round_small(number, ndigits=0):
    mantissa, exponent = frexp10(number)
    mantissa = round(mantissa, ndigits)
    return mantissa*10**exponent