import itertools
import random
import sys
from itertools import combinations
from math import ceil, sqrt, cos, sin, radians, acos
import numpy as np

from pychemia import pcm_log, HAS_PYHULL
from pychemia.utils.mathematics import length_vectors, angle_vectors, wrap2_pmhalf, \
    unit_vector, rotation_matrix_around_axis_angle, angle_vector
from pychemia.utils.mathematics import matrix_from_eig, vector_set_perpendicular
from pychemia.utils.periodic import covalent_radius
from pychemia.core.composition import Composition


class Lattice:
    """
    Routines to create and manipulate the lattice
    The lattice is sufficiently general to account for periodicity in 1, 2 or 3 directions.
    However many routines are only implemented for 3 directions
    The lattice contains 1, 2 or 3 vectors
    """

    def __init__(self, cell=None, periodicity=True):
        """
        Defines the lattice, basically it contains the lattice vectors arranged as rows in the cell array.
        The Lattice class provides methods to compute the reciprocal lattice (That is another Lattice object)
        and  basic routines related to the cell


        >>> cubic = Lattice()
        >>> cubic.lengths # doctest: +SKIP
        array([1., 1., 1.])
        >>> cubic.angles # doctest: +SKIP
        array([90., 90., 90.])

        >>> ortho = Lattice([1, 2, 3])
        >>> ortho.lengths # doctest: +SKIP
        array([1., 2., 3.])
        >>> ortho.angles # doctest: +SKIP
        array([90., 90., 90.])

        >>> bcc = Lattice([[0.5, 0.5, -0.5], [-0.5, 0.5, 0.5], [0.5, -0.5, 0.5]])
        >>> bcc.angles # doctest: +SKIP
        array([109.47122063, 109.47122063, 109.47122063])
        >>> bcc.lengths # doctest: +SKIP
        array([0.8660254, 0.8660254, 0.8660254])

        >>> fcc = Lattice([[0.5, 0.5, 0], [0, 0.5, 0.5], [0.5, 0, 0.5]])
        >>> fcc.lengths # doctest: +SKIP
        array([0.70710678, 0.70710678, 0.70710678])
        >>> fcc.angles # doctest: +SKIP
        array([60., 60., 60.])

        """
        if cell is None:
            cell = np.eye(3)
        else:
            cell = np.array(cell)
        if cell.shape == (3,):
            cell = cell * np.eye(3)
        elif cell.shape == ():
            cell = cell * np.eye(3)

        self._periodicity = None
        self.set_periodicity(periodicity)
        self._dims = sum(self._periodicity)
        assert (np.prod(np.array(cell).shape) == self.periodic_dimensions ** 2)
        self._cell = np.array(cell).reshape((self.periodic_dimensions, self.periodic_dimensions))
        self._lengths = length_vectors(self._cell)
        self._angles = angle_vectors(self._cell, units='deg')
        self._metric = None
        self._inverse = None
        self._limits_for_distance2 = None

    def __str__(self):
        ret = 'Cell='
        for i in range(3):
            for j in range(3):
                ret += "%12.3f " % (self.cell[i, j])
            ret += '\n'
            if i < 2:
                ret += 5 * ' '
        ret += '\n'
        ret += 'Angles: alpha = %12.3f\n' % self.alpha
        ret += '         beta = %12.3f\n' % self.beta
        ret += '        gamma = %12.3f\n' % self.gamma
        ret += '\n'
        ret += 'Lengths:    a = %12.3f\n' % self.a
        ret += '            b = %12.3f\n' % self.b
        ret += '            c = %12.3f\n' % self.c
        return ret

    def cartesian2reduced(self, x):
        return np.dot(x, self.inverse)

    def copy(self):
        """
        Return a copy of the object
        """
        return self.__class__(self._cell, self._periodicity)

    def distance2(self, x1, x2, option='reduced', radius=20, limits=None):

        # Compute the vector from x1 to x2
        dv = np.array(x2) - np.array(x1)

        # If we are not in reduced coordinates,
        # Convert into them
        if option != 'reduced':
            dred = self.cartesian2reduced(dv)
        else:
            dred = dv

        dwrap = wrap2_pmhalf(dred)

        if limits is None:
            limits = np.zeros(3, dtype=int)
            corners = self.get_wigner_seitz_container()

            limits[0] = min(int(ceil(max(1e-14 + abs(np.array([corners[x][0] for x in corners]))))), 5)
            limits[1] = min(int(ceil(max(1e-14 + abs(np.array([corners[x][1] for x in corners]))))), 5)
            limits[2] = min(int(ceil(max(1e-14 + abs(np.array([corners[x][2] for x in corners]))))), 5)

        ret = {}
        for i0 in np.arange(-limits[0], limits[0] + 1):
            for i1 in np.arange(-limits[1], limits[1] + 1):
                for i2 in np.arange(-limits[2], limits[2] + 1):
                    dtot = dwrap + np.array([i0, i1, i2])
                    norm2 = np.dot(np.dot(dtot, self.metric), dtot)
                    if norm2 < radius * radius:
                        ret[(i0, i1, i2)] = {'distance': sqrt(norm2), 'image': dtot}

        return ret

    @staticmethod
    def from_parameters_to_cell(a, b, c, alpha, beta, gamma):
        """
        Create a Lattice using parameters a, b and c as box lengths
        and corresponging angles alpha, beta, gamma

        :param a: (float): *a* lattice length.
        :param b: (float): *b* lattice length.
        :param c: (float): *c* lattice length.
        :param alpha: (float): *alpha* angle in degrees.
        :param beta: (float): *beta* angle in degrees.
        :param gamma: (float): *gamma* angle in degrees.

        :rtype : (Lattice) object
        """

        alpha_radians = radians(alpha)
        beta_radians = radians(beta)
        gamma_radians = radians(gamma)
        val = (cos(alpha_radians) * cos(beta_radians) - cos(gamma_radians)) / (sin(alpha_radians) * sin(beta_radians))
        # Sometimes rounding errors result in values slightly > 1.
        val = val if abs(val) <= 1 else val / abs(val)
        gamma_star = acos(val)
        cell = np.zeros((3, 3))
        cell[0] = [a * sin(beta_radians), 0.0, a * cos(beta_radians)]
        cell[1] = [-b * sin(alpha_radians) * cos(gamma_star),
                   b * sin(alpha_radians) * sin(gamma_star),
                   b * cos(alpha_radians)]
        cell[2] = [0.0, 0.0, float(c)]
        return Lattice(cell)

    def get_brillouin(self):
        return self.reciprocal().get_wigner_seitz()

    def get_path(self):

        assert (self.periodic_dimensions == 3)

        zero3 = np.zeros(3)
        x = self.cell[0, :]
        y = self.cell[1, :]
        z = self.cell[2, :]

        frame = np.array([zero3, x, x + y, y, zero3, z, x + z, x + y + z, y + z, z])

        line1 = np.array([x, x + z])
        line2 = np.array([x + y, x + y + z])
        line3 = np.array([y, y + z])

        return frame, line1, line2, line3

    def get_wigner_seitz(self):
        if HAS_PYHULL:
            from pyhull.voronoi import VoronoiTess
            points = []
            for i, j, k in itertools.product((-1, 0, 1), repeat=3):
                points.append(i * self.cell[0] + j * self.cell[1] + k * self.cell[2])
            tess = VoronoiTess(points)
            ret = []
            for r in tess.ridges:
                if r[0] == 13 or r[1] == 13:
                    ret.append([tess.vertices[i] for i in tess.ridges[r]])
            return ret
        else:
            raise NotImplementedError

    def get_wigner_seitz_container(self):
        """
        Compute the corners of the box that contains the Wigner-Seitz cell

        :return: dict : dictionary with values numpy arrays
        """
        ret = {}
        for i in itertools.product((-1, 1), repeat=3):
            ret[i] = np.dot(self.reciprocal().metric, i * np.diagonal(self.metric))
        return ret

    def minimal_distance(self, x1, x2, option='reduced'):
        distances_dict = self.distance2(x1, x2, option=option)
        mindist = sys.float_info.max
        for k in distances_dict:
            if 0 < distances_dict[k]['distance'] < mindist:
                mindist = distances_dict[k]['distance']
        return mindist

    def minimal_distances(self, reduced1, reduced2):
        """
        Computes a matrix with the minimal distances between
        two sets of points represented as reciprocal coordinates

        :param reduced1: List or array of reduced coordinates for the first
                            set of points
        :param reduced2: Second set of points

        Example:
        >>> import numpy as np
        >>> lattice = Lattice.random_cell('C2')
        >>> r1 = np.random.rand(4, 3)
        >>> r2 = np.random.rand(4, 3)
        >>> dist, close_imgs = lattice.minimal_distances(r1, r2)

        >>> solution = np.zeros((len(r1), len(r2)))

        >>> for i in range(len(r1)):
        ...        for j in range(len(r2)):
        ...                reduced1 = r1[i]
        ...                reduced2 = r2[j]+close_imgs[i, j]
        ...                cartesian1 = lattice.reduced2cartesian(reduced1)
        ...                cartesian2 = lattice.reduced2cartesian(reduced2)
        ...                diff_vector = cartesian2 - cartesian1
        ...                solution[i, j] = np.linalg.norm(diff_vector)

        >>> solution - dist < 1E-5 # doctest: +SKIP
        array([[ True,  True,  True,  True],
               [ True,  True,  True,  True],
               [ True,  True,  True,  True],
               [ True,  True,  True,  True]], dtype=bool)



        """
        # Just in case of one single coordinate
        reduced1, reduced2 = np.atleast_2d(reduced1, reduced2)

        images = np.array([list(i) for i in itertools.product([-1, 0, 1], repeat=3)])

        red2_images = reduced2[:, None, :] + images[None, :, :]

        cartesian1 = self.reduced2cartesian(reduced1)
        cartesian2 = self.reduced2cartesian(red2_images)

        diff_vectors = cartesian2[None, :, :, :] - cartesian1[:, None, None, :]

        distances = np.sum(diff_vectors * diff_vectors, axis=3)

        close_images = np.zeros((len(reduced1), len(reduced2), 3))
        for i in range(len(reduced1)):
            for j in range(len(reduced2)):
                dij = distances[i, j]
                close_images[i, j] = images[dij == min(dij)][0]

        return np.min(distances, axis=2) ** 0.5, close_images

    def distances_in_sphere(self, x1, x2, radius, option='reduced', exclude_out_sphere=True, sort_by_distance=True):
        """
        Returns all the distances between two positions x1 and x2
        taking into account the periodicity of the cell

        :param sort_by_distance:
        :param exclude_out_sphere:
        :param x1:
        :param x2:
        :param radius:
        :param option:
        :return:
        """
        # Compute the vector from x1 to x2
        dv = np.array(x2) - np.array(x1)

        # If the positions are not in reduced coordinates,convert into them
        if option != 'reduced':
            dred = self.cartesian2reduced(dv)
        else:
            dred = dv

        dwrap = wrap2_pmhalf(dred)
        # log.debug('The wrap vector is: %7.3f %7.3f %7.3f' % tuple(dwrap))

        # We need to compute the distances between equivalent faces
        # For that we to compute the perpendicular distance between
        # two faces, the reciprocal lattice produces the inverse of
        # those distances
        recp_len = np.array(self.reciprocal().lengths)
        # The reciprocal (2*pi) is not necessary
        limits = np.ceil(radius * recp_len).astype(int)
        # log.debug('The limits are: %d %d %d' % tuple(limits))

        ret = {}
        # The total size is the product of the number in the 3 directions
        # that is double the limits plus 1 to include the zero
        total_size = np.prod(2 * limits + 1)
        # log.debug('The total size is: %d' % total_size)

        ret['distance'] = np.zeros(total_size)
        ret['image'] = np.zeros((total_size, 3), dtype=int)
        ret['dwrap'] = dwrap
        ret['limits'] = limits

        index = 0
        for i0 in np.arange(-limits[0], limits[0] + 1):
            for i1 in np.arange(-limits[1], limits[1] + 1):
                for i2 in np.arange(-limits[2], limits[2] + 1):
                    dtot = dwrap + np.array([i0, i1, i2])
                    norm2 = np.dot(np.dot(dtot, self.metric), dtot)
                    ret['distance'][index] = sqrt(norm2)
                    ret['image'][index] = np.array([i0, i1, i2])
                    index += 1

        # Exclude distances out of sphere
        if exclude_out_sphere:
            inside_sphere = ret['distance'] <= radius
            if sum(inside_sphere) > 0:
                ret['distance'] = ret['distance'][inside_sphere]
                ret['image'] = ret['image'][inside_sphere]

        if sort_by_distance:
            sorted_indices = np.argsort(ret['distance'])
            ret['distance'] = ret['distance'][sorted_indices]
            ret['image'] = ret['image'][sorted_indices]

        return ret

    @staticmethod
    def random_cell(composition):
        comp = Composition(composition)
        volume = comp.covalent_volume(packing='cubes')

        random.seed()

        # make 3 random lengths
        a = (1.0 + 0.5 * random.random())
        b = (1.0 + 0.5 * random.random())
        c = (1.0 + 0.5 * random.random())

        # now we make 3 random angles
        alpha = 60.0 + 60.0 * random.random()
        beta = 60.0 + 60.0 * random.random()
        gamma = 60.0 + 60.0 * random.random()

        lattice = Lattice().from_parameters_to_cell(a, b, c, alpha, beta, gamma)

        factor = (volume / lattice.volume) ** (1 / 3.0)
        lattice = Lattice().from_parameters_to_cell(factor * a, factor * b, factor * c, alpha, beta, gamma)
        lattice.align_with_axis()
        lattice.align_with_plane()

        return lattice

    def reciprocal(self):
        """
        Return the reciprocal cell

        :rtype : Lattice
        :return:
        """
        return self.__class__(np.linalg.inv(self.cell.T))

    def reduced2cartesian(self, x):
        return np.dot(x, self.cell)

    def set_periodicity(self, periodicity):
        if isinstance(periodicity, bool):
            self._periodicity = 3 * [periodicity]
        elif np.iterable(periodicity):
            assert (len(periodicity) == 3)
            for i in periodicity:
                assert (isinstance(i, bool))
            self._periodicity = list(periodicity)
        else:
            raise ValueError("Periodicity must be a boolean or list instead of: ", periodicity)

    @property
    def cell(self):
        """
        Return the cell as a numpy array

        :return:
        """
        return self._cell

    @property
    def periodic_dimensions(self):
        """
        Return the number of periodic dimensions

        :return: int
        """
        return self._dims

    @property
    def volume(self):
        """
        Computes the volume of the cell (3D),
        area (2D) or generalized volume for
        N dimensions

        :rtype : float
        """
        return abs(np.linalg.det(self.cell))

    @property
    def metric(self):
        if self._metric is None:
            self._metric = np.dot(self.cell, self.cell.T)
        return self._metric

    @property
    def inverse(self):
        if self._inverse is None:
            self._inverse = np.linalg.inv(self.cell)
        return self._inverse

    @property
    def alpha(self):
        return self._angles[(1, 2)]

    @property
    def beta(self):
        return self._angles[(0, 2)]

    @property
    def gamma(self):
        return self._angles[(0, 1)]

    @property
    def angles(self):
        return np.round([self.alpha, self.beta, self.gamma], 14)

    @property
    def a(self):
        return self._lengths[0]

    @property
    def b(self):
        return self._lengths[1]

    @property
    def c(self):
        return self._lengths[2]

    @property
    def lengths(self):
        """
        Return the 3 lenghts of the lattice in the same order as the lattice vectors

        :rtype : list
        """
        return self._lengths

    def align_with_axis(self, axis=0, round_decimals=14):
        a = self.cell[0]
        if axis == 0:
            b = np.array([1, 0, 0])
        elif axis == 1:
            b = np.array([0, 1, 0])
        elif axis == 2:
            b = np.array([0, 0, 1])
        else:
            raise ValueError('Axis must be an integer in (0,1,2)')
        if np.linalg.norm(np.cross(a, b)) < 1E-10:
            return
        c = unit_vector(np.cross(a, b))
        av = angle_vector(a, b)
        rotation_matrix = rotation_matrix_around_axis_angle(c, av)
        self._cell = np.dot(rotation_matrix, self.cell.T).T.round(round_decimals)

    def align_with_plane(self, axis=2, round_decimals=14):
        a = self.cell[1]

        if axis == 0:
            p1 = np.array([0, 1, 0])
            p2 = np.array([0, 0, 1])
        elif axis == 1:
            p1 = np.array([0, 0, 1])
            p2 = np.array([1, 0, 0])
        elif axis == 2:
            p1 = np.array([1, 0, 0])
            p2 = np.array([0, 1, 0])
        else:
            raise ValueError('Axis must be an integer in (0,1,2)')

        if np.linalg.norm(np.cross(p1, a)) < 1E-10:
            return
        c = unit_vector(np.cross(p1, a))
        # print c
        # A vector perpendicular to the plane
        vector_plane = unit_vector(np.cross(p1, p2))
        v1_u = unit_vector(vector_plane)
        v2_u = unit_vector(c)
        proj = np.dot(v1_u, v2_u)
        # print 'Projection', proj
        av = np.arccos(proj)
        # import math
        # print 'Angle=', math.degrees(av)
        rotation_matrix = rotation_matrix_around_axis_angle(p1, -av)
        cell = np.dot(rotation_matrix, self.cell.T).T.round(round_decimals)
        # print '-->',cell[1], vector_plane
        # print 'PROJECTION', np.dot(cell[1], vector_plane)
        if np.abs(np.dot(cell[1], vector_plane)) > 1E-10:
            # print 'Failed projection', np.dot(cell[1], vector_plane)
            # print cell
            rotation_matrix = rotation_matrix_around_axis_angle(p1, av)
            cell = np.dot(rotation_matrix, self.cell.T).T.round(round_decimals)
            if np.dot(cell[1], vector_plane) > 1E-10:
                # print 'Failed projection', np.dot(cell[1], vector_plane)
                # print cell
                pass
        self._cell = cell

    def stretch(self, symbols, rpos, tolerance=1.0, extra=0.1):
        # Dummy array that always will be overwritten
        eigv = np.array([1, 0, 0])
        lattice = self.copy()
        assert len(rpos) == len(symbols)
        natom = len(rpos)
        for i, j in combinations(range(natom), 2):
            ret = lattice.distance2(rpos[i], rpos[j])
            mindist = sys.float_info.max
            for k in ret:
                if 0 < ret[k]['distance'] < mindist:
                    mindist = ret[k]['distance']
                    eigv = ret[k]['image']

            covalent_distance = sum(covalent_radius([symbols[i], symbols[j]]))
            if mindist < tolerance * covalent_distance:
                factor = (tolerance + extra) * covalent_distance / mindist
                v1, v2, v3 = vector_set_perpendicular(eigv)
                matrix_a = matrix_from_eig(v1, v2, v3, factor, 1, 1)
                lattice = Lattice(np.dot(matrix_a, lattice.cell))
        return lattice

    def scale(self, symbols, rpos, tolerance=1.0):
        lattice = self.copy()
        factor = 1.0
        for i in range(len(rpos)):
            # This is to separate each atom from its own image
            covalent_dim = 2.0 * tolerance * covalent_radius(symbols[i])

            this_factor = covalent_dim / lattice.a
            if this_factor > factor:
                factor = this_factor
            this_factor = covalent_dim / lattice.b
            if this_factor > factor:
                factor = this_factor
            this_factor = covalent_dim / lattice.c
            if this_factor > factor:
                factor = this_factor

            for j in range(i + 1, len(rpos)):
                distance = lattice.minimal_distance(rpos[i], rpos[j])
                covalent_dim = tolerance * sum(covalent_radius([symbols[i], symbols[j]]))
                this_factor = covalent_dim / distance
                if this_factor > factor:
                    factor = this_factor
        a = lattice.a
        b = lattice.b
        c = lattice.c
        alpha = lattice.alpha
        beta = lattice.beta
        gamma = lattice.gamma
        return Lattice().from_parameters_to_cell(factor * a, factor * b, factor * c, alpha, beta, gamma)
