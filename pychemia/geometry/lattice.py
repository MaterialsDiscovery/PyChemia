__author__ = 'Guillermo'

import numpy as _np
from math import pi

class Lattice():
    """
    Routines to create and manipulate the lattice
    The lattice is sufficiently general to account for periodicity in 1, 2 or 3 directions.
    However many routines are only implemented for 3 directions
    The lattice contains 1, 2 or 3 vectors
    """

    def __init__(self, cell, periodicity):
        """
        Defines an object lattice that could live
        in 1, 2 or 3 dimensions
        """
        self.cell = _np.array(cell)
        if isinstance(periodicity, bool):
            self.periodicity = 3*[periodicity]
        else:
            self.periodicity = list(periodicity)

    @property
    def periodic_dimensions(self):
        ret = 0
        for i in self.periodicity:
            if i:
                ret += 1
        return ret

    @property
    def volume(self):
        """
        Computes the volume of the cell

        :rtype : float
        """
        return abs(_np.linalg.det(self.cell))

    def reciprocal(self):

        return 2*pi* self.cell.T.I

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
