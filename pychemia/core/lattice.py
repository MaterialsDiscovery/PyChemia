import sys
import random
import itertools
import numpy as _np
from math import ceil, sqrt, cos, sin, radians, acos

from pychemia.utils.mathematics import length_vectors, angle_vectors, wrap2_pmhalf
from composition import Composition
from pychemia import log

__author__ = 'Guillermo Avendano-Franco'


class Lattice():
    """
    Routines to create and manipulate the lattice
    The lattice is sufficiently general to account for periodicity in 1, 2 or 3 directions.
    However many routines are only implemented for 3 directions
    The lattice contains 1, 2 or 3 vectors
    """

    def __init__(self, cell=None, periodicity=True):
        """
        Defines an object lattice that could live
        in arbitrary dimensions
        """
        if cell is None:
            cell = _np.eye(3)
        self._periodicity = None
        self.set_periodicity(periodicity)
        self._dims = sum(self._periodicity)
        assert (_np.prod(_np.array(cell).shape) == self.periodic_dimensions ** 2)
        self._cell = _np.array(cell).reshape((self.periodic_dimensions, self.periodic_dimensions))
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
        return _np.dot(x, self.inverse)

    def copy(self):
        """
        Return a copy of the object
        """
        return self.__class__(self._cell, self._periodicity)

    def distance2(self, x1, x2, option='reduced', radius=20, limits=None):

        # Compute the vector from x1 to x2
        dv = _np.array(x2) - _np.array(x1)

        # If we are not in reduced coordinates,
        # Convert into them
        if option != 'reduced':
            dred = self.cartesian2reduced(dv)
        else:
            dred = dv

        dwrap = wrap2_pmhalf(dred)

        if limits is None:
            limits = _np.zeros(3, dtype=int)
            corners = self.get_wigner_seitz_container()

            limits[0] = min(int(ceil(max(1e-14 + abs(_np.array([corners[x][0] for x in corners]))))), 5)
            limits[1] = min(int(ceil(max(1e-14 + abs(_np.array([corners[x][1] for x in corners]))))), 5)
            limits[2] = min(int(ceil(max(1e-14 + abs(_np.array([corners[x][2] for x in corners]))))), 5)

        ret = {}
        for i0 in _np.arange(-limits[0], limits[0] + 1):
            for i1 in _np.arange(-limits[1], limits[1] + 1):
                for i2 in _np.arange(-limits[2], limits[2] + 1):
                    dtot = dwrap + _np.array([i0, i1, i2])
                    norm2 = _np.dot(_np.dot(dtot, self.metric), dtot)
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
        cell = _np.zeros((3, 3))
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

        zero3 = _np.zeros(3)
        x = self.cell[0, :]
        y = self.cell[1, :]
        z = self.cell[2, :]

        frame = _np.array([zero3, x, x + y, y, zero3, z, x + z, x + y + z, y + z, z])

        line1 = _np.array([x, x + z])
        line2 = _np.array([x + y, x + y + z])
        line3 = _np.array([y, y + z])

        return frame, line1, line2, line3

    def get_wigner_seitz(self):

        from pyhull.voronoi import VoronoiTess
        import itertools

        points = []
        for i, j, k in itertools.product((-1, 0, 1), repeat=3):
            points.append(i * self.cell[0] + j * self.cell[1] + k * self.cell[2])
        tess = VoronoiTess(points)
        ret = []
        for r in tess.ridges:
            if r[0] == 13 or r[1] == 13:
                ret.append([tess.vertices[i] for i in tess.ridges[r]])
        return ret

    def get_wigner_seitz_container(self):
        """
        Compute the corners of the box that contains the Wigner-Seitz cell

        :return: dict : dictionary with values numpy arrays
        """
        ret = {}
        for i in itertools.product((-1, 1), repeat=3):
            ret[i] = _np.dot(self.reciprocal().metric, i * _np.diagonal(self.metric))
        return ret

    def minimal_distance(self, x1, x2, option='reduced'):
        distances_dict = self.distance2(x1, x2, option=option)
        mindist = sys.float_info.max
        for k in distances_dict:
            if 0 < distances_dict[k]['distance'] < mindist:
                mindist = distances_dict[k]['distance']
        return mindist

    def minimal_distances(self, red_coords1, red_coords2):
        """
        Computes a matrix with the minimal distances between
        two sets of points represented as reciprocal coordinates

        :param red_coords1: List or array of reduced coordinates for the first
                            set of points
        :param red_coords2: Second set of points
        """
        # Just in case of one single coordinate
        red_coords1, red_coords2 = _np.atleast_2d(red_coords1, red_coords2)

        images = _np.array([list(i) for i in itertools.product([-1, 0, 1], repeat=3)])

        red2_images = red_coords2[:, None, :] + images[None, :, :]

        cart_coords1 = self.reduced2cartesian(red_coords1)
        cart_coords2 = self.reduced2cartesian(red2_images)

        diff_vectors = cart_coords2[None, :, :, :] - cart_coords1[:, None, None, :]
        return _np.min(_np.sum(diff_vectors ** 2, axis=3), axis=2) ** 0.5

    def distances_in_sphere(self, x1, x2, radius, option='reduced', exclude_out_sphere=True, sort_by_distance=True):
        """
        Returns all the distances between two positions x1 and x2
        taking into account the periodicity of the cell

        :param x1:
        :param x2:
        :param radius:
        :param option:
        :return:
        """
        # Compute the vector from x1 to x2
        dv = _np.array(x2) - _np.array(x1)

        # If the positions are not in reduced coordinates,convert into them
        if option != 'reduced':
            dred = self.cartesian2reduced(dv)
        else:
            dred = dv

        dwrap = wrap2_pmhalf(dred)
        log.debug('The wrap vector is: %7.3f %7.3f %7.3f' % tuple(dwrap))

        # We need to compute the distances between equivalent faces
        # For that we to compute the perpendicular distance between
        # two faces, the reciprocal lattice produces the inverse of
        # those distances
        recp_len = _np.array(self.reciprocal().lengths)
        # The reciprocal (2*pi) is not necessary
        limits = _np.ceil(radius * recp_len).astype(int)
        log.debug('The limits are: %d %d %d' % tuple(limits))

        ret = {}
        # The total size is the product of the number in the 3 directions
        # that is double the limits plus 1 to include the zero
        total_size = _np.prod(2 * limits + 1)
        log.debug('The total size is: %d' % total_size)

        ret['distance'] = _np.zeros(total_size)
        ret['image'] = _np.zeros((total_size, 3), dtype=int)
        ret['dwrap'] = dwrap
        ret['limits'] = limits

        index = 0
        for i0 in _np.arange(-limits[0], limits[0] + 1):
            for i1 in _np.arange(-limits[1], limits[1] + 1):
                for i2 in _np.arange(-limits[2], limits[2] + 1):
                    dtot = dwrap + _np.array([i0, i1, i2])
                    norm2 = _np.dot(_np.dot(dtot, self.metric), dtot)
                    ret['distance'][index] = sqrt(norm2)
                    ret['image'][index] = _np.array([i0, i1, i2])
                    index += 1

        # Exclude distances out of sphere
        if exclude_out_sphere:
            inside_sphere = ret['distance'] <= radius
            ret['distance'] = ret['distance'][inside_sphere]
            ret['image'] = ret['image'][inside_sphere]

        if sort_by_distance:
            sorted_indices = _np.argsort(ret['distance'])
            ret['distance'] = ret['distance'][sorted_indices]
            ret['image'] = ret['image'][sorted_indices]

        return ret

    def plot(self, points=None):
        from mayavi import mlab

        tube_rad = max(self.lengths) / 100.0

        frame, line1, line2, line3 = self.get_path()
        for i, j, k in [[frame[:, 0], frame[:, 1], frame[:, 2]],
                        [line1[:, 0], line1[:, 1], line1[:, 2]],
                        [line2[:, 0], line2[:, 1], line2[:, 2]],
                        [line3[:, 0], line3[:, 1], line3[:, 2]]]:
            mlab.plot3d(i, j, k, tube_radius=tube_rad, color=(1, 1, 1), tube_sides=24, transparent=True, opacity=0.5)

        if points is not None:
            ip = _np.array(points)
            mlab.points3d(ip[:, 0], ip[:, 1], ip[:, 2], tube_rad * _np.ones(len(ip)), scale_factor=1)

        return mlab.pipeline

    def plot_wigner_seitz(self, scale=1):
        import itertools
        from tvtk.api import tvtk

        ws = _np.array(self.get_wigner_seitz())
        points = _np.array([])
        for i in ws:
            for j in i:
                points = _np.concatenate((points, j))
        points = scale * (points.reshape(-1, 3))

        index = 0
        # triangles = index + _np.array(list(itertools.combinations(range(len(ws[0])), 3)))
        #scalars = _np.ones(len(ws[0]))
        #for i in ws[1:]:
        #    index += len(i)
        #    triangles = _np.concatenate((triangles, index + _np.array(list(itertools.combinations(range(len(i)), 3)))))
        #    scalars = _np.concatenate((scalars, _np.random.random() * _np.ones(len(i))))
        #scalars = _np.ones(len(ws[0]))
        scalars = None
        triangles = None
        for i in ws:
            print i
            iscalars = _np.ones(len(i))
            itriangles = index + _np.array(list(itertools.combinations(range(len(i)), 3)))
            print iscalars
            print itriangles
            if triangles is None:
                triangles = itriangles
            else:
                triangles = _np.concatenate((triangles, itriangles))
            if scalars is None:
                scalars = iscalars
            else:
                scalars = _np.concatenate((scalars, iscalars))
            index += len(i)

        #print triangles
        #print scalars
        # The TVTK dataset.
        mesh = tvtk.PolyData(points=points, polys=triangles)
        mesh.point_data.scalars = scalars
        mesh.point_data.scalars.name = 'scalars'

        pipeline = self.plot()
        pipeline.surface(mesh, color=(0.9, 0.1, 0.1), opacity=0.2)
        return pipeline

    @staticmethod
    def random_cell(composition):

        if isinstance(composition, dict):
            comp = Composition(composition)
        elif isinstance(composition, Composition):
            comp = composition
        else:
            raise ValueError('Wrong composition value')

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

        return lattice

    def reciprocal(self):
        """
        Return the reciprocal cell

        :rtype : Lattice
        :return:
        """
        return self.__class__(_np.linalg.inv(self.cell.T))

    def reduced2cartesian(self, x):
        return _np.dot(x, self.cell)

    def set_periodicity(self, periodicity):
        if isinstance(periodicity, bool):
            self._periodicity = 3 * [periodicity]
        elif _np.iterable(periodicity):
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
        return abs(_np.linalg.det(self.cell))

    @property
    def metric(self):
        if self._metric is None:
            self._metric = _np.dot(self.cell, self.cell.T)
        return self._metric

    @property
    def inverse(self):
        if self._inverse is None:
            self._inverse = _np.linalg.inv(self.cell)
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
        return self.alpha, self.beta, self.gamma

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

