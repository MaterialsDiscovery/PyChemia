import numpy as _np
from math import pi
from pychemia.utils.mathematics import length_vectors, angle_vectors

__author__ = 'Guillermo Avendano-Franco'


class Lattice():
    """
    Routines to create and manipulate the lattice
    The lattice is sufficiently general to account for periodicity in 1, 2 or 3 directions.
    However many routines are only implemented for 3 directions
    The lattice contains 1, 2 or 3 vectors
    """

    def __init__(self, cell, periodicity=True):
        """
        Defines an object lattice that could live
        in arbitrary dimensions
        """
        if isinstance(periodicity, bool):
            self._periodicity = 3*[periodicity]
        elif isinstance(periodicity, list):
            self._periodicity = list(periodicity)
        else:
            raise ValueError("periodicity must be a boolean or list")

        self._dims = sum(self._periodicity)
        assert(_np.prod(_np.array(cell).shape) == self.periodic_dimensions**2)
        self._cell = _np.array(cell).reshape((self.periodic_dimensions, self.periodic_dimensions))
        self._lengths = length_vectors(self._cell)
        self._angles = angle_vectors(self._cell)
        self._metric = None
        self._inverse = None

    @staticmethod
    def parameters2cell(a, b, c, alpha, beta, gamma):
        pass

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
            self._inverse = self.cell.I
        return self._inverse

    @property
    def reciprocal(self):
        """
        Return the reciprocal cell

        :rtype : Lattice
        :return:
        """
        return self.__class__(2*pi*self.cell.T.I)

    def copy(self):
        """
        Return a copy of the object
        """
        return self.__class__(self._cell, self._periodicity)

    def get_path(self):

        assert(self.periodic_dimensions == 3)

        cero = _np.zeros(3)
        x = self.cell[0, :]
        y = self.cell[1, :]
        z = self.cell[2, :]

        frame = _np.array([cero, x, x+y, y, cero, z, x+z, x+y+z, y+z, z])

        line1 = _np.array([x, x+z])
        line2 = _np.array([x+y, x+y+z])
        line3 = _np.array([y, y+z])

        return frame, line1, line2, line3

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
        return self._lengths

