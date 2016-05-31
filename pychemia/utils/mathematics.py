from __future__ import print_function
import itertools
import math
from collections import OrderedDict
from fractions import gcd
from math import cos, sin, sqrt
import numpy as np


def length_vector(v):
    """
    Returns the length of a vector 'v' in arbitrary number of dimensions

    :param v: (list, numpy.ndarray) Vector to compute length

    :rtype : (float) The lenght of the vector

    Examples:
    >>> length_vector([1, 2, 3])
    3.7416573867739413

    """
    return np.linalg.norm(v)


def length_vectors(m):
    """
    Returns the lengths of several vectors
    arranged as rows in a MxN matrix

    :param m: numpy.ndarray

    :rtype : numpy.ndarray

    Examples:
    >>> length_vectors([[1, 2, 3], [4, 5, 6], [7, 8, 9], [1, 0, 0], [0, 0, 2]])
    array([  3.74165739,   8.77496439,  13.92838828,   1.        ,   2.        ])

    """
    m = np.array(m)
    return np.linalg.norm(m, axis=1)


def unit_vector(v):
    """
    Returns the unit vector of the vector.
    Arbitrary number of dimensions

    :param v: list, numpy.array

    :rtype : numpy.ndarray

    Examples:
    >>> a = unit_vector([1, 2, 3])
    >>> a
    array([ 0.26726124,  0.53452248,  0.80178373])

    >>> length_vector(a)
    1.0

    """
    if length_vector(np.array(v, dtype=float)) < 1E-10:
        raise ValueError('Vector is null')
    return np.array(v) / length_vector(np.array(v, dtype=float))


def unit_vectors(m):
    """
    Returns the unit vectors of a set
    of vectors arranged as rows in MxN matrix

    :param m: numpy.ndarray

    :rtype : numpy.ndarray

    Example:
    >>> from pychemia.utils.mathematics import *
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
    m = np.array(m)
    return m / (np.linalg.norm(m, axis=1)[:, np.newaxis])


def angle_vector(v1, v2, units='rad'):
    """
    Returns the angle in radians (default) or degrees
    between vectors 'v1' and 'v2'::

    :param v1: (list, numpy.ndarray)
    :param v2: (list, numpy.ndarray)
    :param units: (str) : 'rad' (default) Radians, 'deg' Degrees

    :rtype : float

    Examples:
    >>> angle_vector([1, 0, 0], [0, 1, 0])
    1.5707963267948966
    >>> angle_vector([1, 0, 0], [1, 0, 0])
    0.0
    >>> angle_vector([1, 0, 0], [-1, 0, 0])
    3.1415926535897931
    >>> angle_vector([1, 0, 0], [0, 1, 0], units='deg')
    90.0
    >>> angle_vector([1, 0, 0], [-1, 0, 0], units='deg')
    180.0

    """
    assert (units in ['rad', 'deg'])

    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    angle = np.arccos(np.dot(v1_u, v2_u))
    if np.isnan(angle):
        if (v1_u == v2_u).all():
            return 0.0
        else:
            return np.pi
    if units == 'rad':
        return angle
    elif units == 'deg':
        return 180.0 * angle / np.pi


def angle_between_vectors(a, b):
    assert a.shape == b.shape

    norms = (np.linalg.norm(a, axis=1) * np.linalg.norm(b, axis=1))
    norms[norms == 0] = 1
    dots = np.sum(a * b, axis=1) / norms
    dots = np.round(dots, 15)
    dots[dots > 1] = 1.0
    norms = (np.linalg.norm(a, axis=1) * np.linalg.norm(b, axis=1))
    dots[norms == 0] = 1.0
    return np.arccos(dots)


def angle_vectors(m, units='rad'):
    """
    Returns all the angles for all the
    vectors arranged as rows in matrix 'm'

    :param m: (numpy.ndarray)
    :param units: (str) : 'rad' Radians, 'deg' Degrees

    :rtype : numpy.ndarray

    Examples:
    >>> import pprint
    >>> a = angle_vectors([[1, 2, 3], [4, 5, 6], [7, 8, 9], [1, 0, 0], [0, 0, 2]])
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
    for i in itertools.combinations(range(len(m)), 2):
        ret[i] = angle_vector(m[i[0]], m[i[1]], units=units)
    return ret


def distance(v1, v2):
    """
    Return the vector v2-v1, the vector going from v1 to v2
    and the magnitude of that vector.

    :param v1: (list, numpy.ndarray)
    :param v2: (list, numpy.ndarray)

    :rtype : tuple

    Examples:
    >>> distance([0, 0, 0, 1], [1, 0, 0, 0])
    (array([ 1,  0,  0, -1]), 1.4142135623730951)
    >>> distance([-1, 0, 0], [1, 0, 0])
    (array([2, 0, 0]), 2.0)

    """
    ret = np.array(v2) - np.array(v1)
    return ret, length_vector(ret)


def distances(m):
    """
    Return all the distances for all possible combinations of the row vectors in matrix m

    :param m: (list, numpy.ndarray)

    :rtype : dict

    Example:
    >>> import pprint
    >>> pprint.pprint(distances([[1, 2, 3], [4, 5, 6], [7, 8, 9], [1, 0, 0], [0, 0, 2]]))
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
    for i in itertools.combinations(range(len(m)), 2):
        ret[i] = distance(m[i[0]], m[i[1]])
    return ret


def wrap2_pmhalf(x):
    """
    Wraps a number or array in the interval ]-1/2, 1/2] values = -1/2 will be wrapped  to 1/2

    :param x: (float) The number to be wrapped in the interval (-1/2, 1/2]

    Examples:
    >>> wrap2_pmhalf(-0.5)
    0.5
    >>> wrap2_pmhalf(0.0)
    0.0
    >>> wrap2_pmhalf([-0.75, -0.5, -0.25, 0.0, 0.25, 0.5, 0.75])
    array([ 0.25,  0.5 , -0.25,  0.  ,  0.25,  0.5 , -0.25])
    >>> wrap2_pmhalf([[-0.75, -0.5, -0.25], [0.25, 0.5, 0.75]])
    array([[ 0.25,  0.5 , -0.25],
           [ 0.25,  0.5 , -0.25]])

    """

    def wrap(num):
        tol12 = 1e-12
        if num > 0:
            ret = (num + 0.5 - tol12) % 1.0 - 0.5 + tol12
        else:
            ret = -(-(num - 0.5 - tol12) % 1.0) + 0.5 + tol12
        for y in [-0.25, 0.0, 0.25, 0.5]:
            ret = (lambda num2: y if abs(y - num2) < tol12 else num2)(ret)
        return ret

    if np.iterable(x):
        vec = np.vectorize(wrap)
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
        other = unit_vector(np.random.random_sample(3))
        if np.abs(np.dot(v1, other)) > 0.05:
            v2 = unit_vector(np.cross(v1, other))
            v3 = unit_vector(np.cross(v1, v2))
            break
        else:
            continue
    # print _np.dot(v1, v2)
    # print _np.dot(v1, v3)
    # print _np.dot(v2, v3)

    # assert (_np.abs(_np.dot(v1, v2)) < 1E-15)
    # assert (_np.abs(_np.dot(v1, v3)) < 1E-15)
    # assert (_np.abs(_np.dot(v2, v3)) < 1E-15)
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
    matrixp = np.vstack((v1, v2, v3)).T
    matrixd = np.diag([lam1, lam2, lam3])
    matrixpinv = np.linalg.inv(matrixp)
    matrixa = np.dot(matrixp, np.dot(matrixd, matrixpinv))
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

    # Integral from -\infty to a
    val_floor = 0.5 * (1 + math.erf((a - mu) / (sigma * math.sqrt(2.0))))

    # Integral from -\infty to b
    val_ceil = 0.5 * (1 + math.erf((b - mu) / (sigma * math.sqrt(2.0))))

    return val_ceil - val_floor


def frexp10(x):
    exp = int(math.floor(math.log10(abs(x))))
    return x / 10 ** exp, exp


def round_small(number, ndigits=0):
    mantissa, exponent = frexp10(number)
    mantissa = round(mantissa, ndigits)
    return mantissa * 10 ** exponent


def sieve_atkin(limit):
    ret = [2, 3]
    sieve = [False] * (limit + 1)
    for x in range(1, int(math.sqrt(limit)) + 1):
        for y in range(1, int(math.sqrt(limit)) + 1):
            n = 4 * x * x + y * y
            if n <= limit and (n % 12 == 1 or n % 12 == 5):
                sieve[n] = not sieve[n]
            n = 3 * x * x + y * y
            if n <= limit and n % 12 == 7:
                sieve[n] = not sieve[n]
            n = 3 * x * x - y * y
            if x > y and n <= limit and n % 12 == 11:
                sieve[n] = not sieve[n]
    for x in range(5, int(math.sqrt(limit))):
        if sieve[x]:
            for y in range(x * x, limit + 1, x * x):
                sieve[y] = False
    for p in range(5, limit):
        if sieve[p]:
            ret.append(p)
    return ret


def trial_division(n):
    """
    Return a list of the prime factors for a natural number
    uses the Sieve of Atkin as a list of primes

    :param n: (int) A natural number

    :rtype : (list)
    """
    if n == 1:
        return [1]
    primes = sieve_atkin(int(n ** 0.5) + 1)
    prime_factors = []

    for p in primes:
        if p * p > n:
            break
        while n % p == 0:
            prime_factors.append(p)
            n /= p
    if n > 1:
        prime_factors.append(n)

    return prime_factors


def lcm(a, b):
    """
    Return the Least Common Multiple
    the smallest positive integer that is divisible by both a and b
    :param a: (int)
    :param b: (int)
    :return: (int)

    :rtype : (int)
    """
    return a * b / gcd(a, b)


def shortest_triple_set(n):
    """
    Return the smallest three numbers (a,b,c) such as
    a*b*c=n
    And a+b+c is as small as possible

    :param n:
    :return:
    """
    # First Compute the prime factors
    prime_factors = trial_division(n)

    if len(prime_factors) == 3:
        # No choice return the three numbers
        return prime_factors
    elif len(prime_factors) < 3:
        # If there are less than 3 complete the set with ones
        while len(prime_factors) % 3 != 0:
            prime_factors = [1] + prime_factors
        return prime_factors
    else:
        factors = np.array(prime_factors)
        while len(factors) > 3:
            # print factors
            # Complete a multiple of 6 and sum folding lowest with highest
            while len(factors) % 6 != 0:
                factors = np.concatenate(([1], factors))
            # Take the first half
            low = factors[:int(len(factors) / 2)]
            # take the second half and invert the order
            high = factors[int(len(factors) / 2):][::-1]
            # Sum both arrays and sort them before reenter
            factors = np.sort(low * high)
        return [int(x) for x in factors]


def rotation_matrix_around_axis_angle(axis, theta):
    """
    Return the rotation matrix needed to rotate any vector
    around the 'axis' an angle of 'theta' radians
    """
    # Given a unit vector u = (ux, uy, uz), where ux**2 + uy**2 + uz**2 = 1,
    u = unit_vector(axis)
    # with ux is the cross product matrix
    ux = np.array([[0, -u[2], u[1]], [u[2], 0, -u[0]], [-u[1], u[0], 0]])
    # uxu is the tensor product of u
    uxu = np.tensordot(u, u.T, axes=0)
    # This is a matrix form of Rodrigues 'rotation formula'
    return np.cos(theta) * np.identity(3) + np.sin(theta) * ux + (1 - np.cos(theta)) * uxu


def rotate_towards_axis(vector, axis, theta=None, fraction=None):
    # If the vector is already parallel to the axis, do nothing
    if np.abs(angle_vector(vector, axis)) < 1E-10 or np.abs(angle_vector(vector, axis) - np.pi) < 1E-10:
        return vector

    # Create a unitary vector perpendicular to the plane created by vector and axis
    uv = unit_vector(np.cross(vector, axis))

    # Use the rodrigues formula to rotate around the vector perpendicular
    if theta is not None:
        m = rotation_matrix_around_axis_angle(uv, theta)
        return np.dot(m, vector)

    if fraction is not None:
        theta = angle_vector(vector, axis)
        m = rotation_matrix_around_axis_angle(uv, fraction * theta)
        return np.dot(m, vector)


def rotation_x(theta):
    """
    Create a rotation matrix around the 'x' axis

    :param theta: (float) Angle in radians
    :return: (numpy.ndarray)

    Examples:
    >>> import numpy as np
    >>> m = rotation_x(np.pi/3)
    >>> np.all(np.round(np.dot(m.T,m),15)==np.eye(3))
    True

    """
    return np.array([[1, 0, 0],
                     [0, np.cos(theta), -np.sin(theta)],
                     [0, np.sin(theta), np.cos(theta)]])


def rotation_y(theta):
    """
    Create a rotation matrix around the 'y' axis

    :param theta: (float) Angle in radians
    :return: (numpy.ndarray)

    Example:
    >>> import numpy as np
    >>> m = rotation_y(np.pi/3)
    >>> np.all(np.round(np.dot(m.T, m), 15)==np.eye(3))
    True

    """
    return np.array([[np.cos(theta), 0, np.sin(theta)],
                     [0, 1, 0],
                     [-np.sin(theta), 0, np.cos(theta)]])


def rotation_z(theta):
    """
    Create a rotation matrix around the 'z' axis

    :param theta: (float) Angle in radians
    :return: (numpy.ndarray)

    Examples:
    >>> import numpy as np
    >>> m = rotation_z(np.pi/3)
    >>> np.all(np.round(np.dot(m.T,m),15)==np.eye(3))
    True

    """
    return np.array([[np.cos(theta), -np.sin(theta), 0],
                     [np.sin(theta), np.cos(theta), 0],
                     [0, 0, 1]])


def apply_rotation(vector, theta_x, theta_y, theta_z):
    """
    Apply a rotation matrix to a vector by succesive rotations around
    the three axis 'x', 'y' and 'z'

    :param vector:
    :param theta_x: (float) Angle in radians
    :param theta_y: (float) Angle in radians
    :param theta_z: (float) Angle in radians
    :return: (numpy.ndarray)

  Example:
    >>> a=apply_rotation([0.1, 0.2, 0.3], 3.1415/3, 3.1415/4, 3.1415/5)
    >>> b=apply_rotation(a, -3.1415/3, 0, 0)
    >>> c=apply_rotation(b, 0, -3.1415/4, 0)
    >>> d=apply_rotation(c, 0, 0, -3.1415/5)
    >>> d
    array([ 0.1,  0.2,  0.3])

    """
    return np.round(np.dot(rotation_x(theta_x), np.dot(rotation_y(theta_y), np.dot(rotation_z(theta_z), vector))), 14)


def rotation_matrix_weave(axis, theta, mat=None):
    try:
        from scipy import weave
    except ImportError:
        raise NotImplementedError
    if mat is None:
        mat = np.eye(3, 3)

    support = "#include <math.h>"
    code = """
        double x = sqrt(axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2]);
        double a = cos(theta / 2.0);
        double b = -(axis[0] / x) * sin(theta / 2.0);
        double c = -(axis[1] / x) * sin(theta / 2.0);
        double d = -(axis[2] / x) * sin(theta / 2.0);

        mat[0] = a*a + b*b - c*c - d*d;
        mat[1] = 2 * (b*c - a*d);
        mat[2] = 2 * (b*d + a*c);

        mat[3*1 + 0] = 2*(b*c+a*d);
        mat[3*1 + 1] = a*a+c*c-b*b-d*d;
        mat[3*1 + 2] = 2*(c*d-a*b);

        mat[3*2 + 0] = 2*(b*d-a*c);
        mat[3*2 + 1] = 2*(c*d+a*b);
        mat[3*2 + 2] = a*a+d*d-b*b-c*c;
    """

    weave.inline(code, ['axis', 'theta', 'mat'], support_code=support, libraries=['m'])

    return mat


def rotation_matrix_numpy(axis, theta):
    mat = np.eye(3, 3)
    axis /= sqrt(np.dot(axis, axis))
    a = cos(theta / 2.)
    b, c, d = -axis * sin(theta / 2.)

    return np.array([[a * a + b * b - c * c - d * d, 2 * (b * c - a * d), 2 * (b * d + a * c)],
                     [2 * (b * c + a * d), a * a + c * c - b * b - d * d, 2 * (c * d - a * b)],
                     [2 * (b * d - a * c), 2 * (c * d + a * b), a * a + d * d - b * b - c * c]])


def projector(u, v):
    """
    Computes the projector of 'v' over 'u'
    A vector in the direction of 'u' with a
    magnitude of the projection of 'v' over 'u'

    :param u:
    :param v:
    :return:

    Example:
    >>> projector([0.1, 0.2, 0.3], [0.3, 0.2, 0.1])
    array([ 0.07142857,  0.14285714,  0.21428571])
    >>> projector([1, 0, 0], [0, 2, 0])
    array([ 0.,  0.,  0.])

    """
    u = np.array(u, dtype=float)
    v = np.array(v, dtype=float)
    return np.dot(v, u) / np.dot(u, u) * np.array(u)


def gram_smith(m):
    """
    Create a Gram-Smith ortoganalized
    matrix, using the first vector of
    matrix 'm' to build the ortogonal
    set of vectors

    :param m:
    :return:

    Example:
    >>> import numpy as np
    >>> o = gram_smith(np.random.rand(3,3))
    >>> np.round(np.abs(np.linalg.det(o)),10)==1.0
    True

    """
    m = np.array(m)
    ret = np.zeros((len(m[0]), len(m[0])))
    ret[0] = m[0] / np.linalg.norm(m[0])
    for k in range(1, len(m[0])):
        ret[k] = m[k]
        for j in range(0, k):
            ret[k] -= projector(ret[j], m[k])
        ret[k] /= np.linalg.norm(ret[k])
    return ret


def gram_smith_qr(ndim=3):
    """
    Create a Gram-Smith orthoganalized matrix.
    The argument is the dimension of the matrix
    and uses a random matrix and QR decomposition
    to build the orthogonal matrix

    :param ndim: Dimension of the matrix
    :return:

    Example:
    >>> import numpy as np
    >>> o = gram_smith_qr(3)
    >>> np.round(np.abs(np.linalg.det(o)),10)==1.0
    True

    """
    matrix_a = np.random.rand(ndim, ndim)
    while np.linalg.det(matrix_a) < 1E-5:
        matrix_a = np.random.rand(ndim, ndim)
    return np.linalg.qr(matrix_a)[0]


def cartesian_to_spherical(xyz):
    """
    Convert a numpy array (Nx3) into
    a numpy array (Nx3) where the first
    column is the magnitude of the vector
    and second and third are spherical angles

    We use the mathematical notation
    (radial, azimuthal, polar)

    :param xyz: numpy.array
    :return: numpy.array
    """
    xyz = np.array(xyz).reshape((-1, 3))
    spherical = np.zeros(xyz.shape)
    xy = xyz[:, 0] ** 2 + xyz[:, 1] ** 2
    spherical[:, 0] = np.sqrt(xy + xyz[:, 2] ** 2)
    spherical[:, 1] = np.arctan2(xyz[:, 1], xyz[:, 0])
    spherical[:, 2] = np.arctan2(np.sqrt(xy), xyz[:, 2])  # for elevation angle defined from Z-axis down
    # spherical[:, 2] = np.arctan2(xyz[:, 2], np.sqrt(xy))   # for elevation angle defined from XY-plane up
    return spherical


def spherical_to_cartesian(spherical):
    """
    Convert a numpy array (Nx3) from
    spherical coordinates to cartesian
    coordinates

    We use the mathematical notation
    (radial, azimuthal, polar)

    :param spherical: numpy.array
    :return: numpy.array
    """
    spherical = np.array(spherical).reshape((-1, 3))
    xyz = np.zeros(spherical.shape)

    xyz[:, 0] = spherical[:, 0] * np.cos(spherical[:, 1]) * np.sin(spherical[:, 2])
    xyz[:, 1] = spherical[:, 0] * np.sin(spherical[:, 1]) * np.sin(spherical[:, 2])
    xyz[:, 2] = spherical[:, 0] * np.cos(spherical[:, 2])
    return xyz


def rotation_ndim(n, t, indices):
    rot = np.eye(n)
    rot[indices[0], indices[0]] = np.cos(t)
    rot[indices[1], indices[1]] = np.cos(t)
    rot[indices[0], indices[1]] = -np.sin(t)
    rot[indices[1], indices[0]] = np.sin(t)
    return rot


def generalized_euler_angles(m):
    # Make a copy of the original matrix
    mp = np.array(m)
    # The dimension of the matrix
    n = m.shape[0]
    # The angles is an ordered list of operations
    # Rotations form a non-abelian group
    angles = OrderedDict()
    for i in itertools.combinations(range(n), 2):
        # Compute the angle to align the plane i
        # into the canonical base
        theta = np.arctan(mp[i[0], i[1]] / mp[i[0], i[0]])
        angles[i] = theta
        # Compute a rotation matrix for plane i
        rot = rotation_ndim(n, theta, i)
        # Apply the rotation left side
        mp = np.array(np.dot(rot, mp))
        retm = rotation_ndim(n, angles)

    # angles is a Ordered dictionary of Generalized Euler Angles
    # mp is a matrix that should be close to Identity
    # retm is the return matrix, the matrix that should be
    # close to the original one
    return angles, mp, retm


def rotation_matrix_ndim(n, angles):
    ret = np.eye(n)
    for i in angles.keys()[::-1]:
        rot = rotation_ndim(n, angles[i], i)
        ret = np.dot(rot.T, ret)
    return ret
