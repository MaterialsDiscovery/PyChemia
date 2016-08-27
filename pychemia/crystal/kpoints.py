from math import sqrt, ceil
import numpy as np
from pychemia.utils.serializer import PyChemiaJsonable, generic_serializer

"""
Definition of the set of k-points in reciprocal space
"""

__author__ = "Guillermo Avendano-Franco"
__copyright__ = "Copyright 2014"
__version__ = "1.1"
__maintainer__ = "Guillermo Avendano-Franco"
__email__ = "gtux.gaf@gmail.com"
__status__ = "Development"
__date__ = "November 5, 2014"


class KPoints(PyChemiaJsonable):
    """
    Defines an object that contains information about kpoints mesh
    There are two mutually exclusive options:

    Explicit array of kpoints with associated weights OR
    a specification of a grid with an optional set of shifts.

    The default is one single kpoint in gamma
    """

    def __init__(self, kmode='gamma',
                 kpoints_list=None, weights=None,
                 shifts=None, grid=None,
                 kvertices=None, intermediates=None):
        """
        The KPoints class is designed to  store the description of k-points in a code agnostic way.

        There are 5 possible kmodes to describe a set of k-points, currently those modes are not interchangable,
        but that could change in the future.

        The kmodes are 'cartesian', 'reduced', 'gamma', 'monkhorst-pack' and 'path'.

        Kmode 'cartesian' and 'reduced' uses explicit declaration of a grid of k-points.
        In the case of 'cartesian' the points are vectors in k-space.
        In the case of 'reduced' the points are the factors for the vectors defining the reciprocal lattice.
        The variable 'kpoints_list' is a numpy array with a shape of ('nkpt', 3)
        The variable 'weigths' is a numpy array with the weights associated to each point.
        by default all the weights are equal to  1.

        Kmode 'gamma' and 'monkhorst-pack' uses an uniform grid to describe the k-points mesh.
        The variable 'grid' describes the number of points on each dimension an it is a numpy array
        with 3 values.
        The variable 'shifts' describe the shift done to the grid according to the mode selected.
        The internal value 'nkpt' will be simply the product of the three values of the grid.
        Typically ab-initio packages will try to exclude from that grid equivalent k-points due to
        symmetry considerations.

        Kmode 'path' is specially useful for band-structure calculations, in this case the variable
        'kvertices' contains all the vertices, usually high symmetry points and the variable
        'intermediates' contains the number of k-points created between two consecutive vertices.

        When a Kpoints object is created without any arguments, the default will be a kmode='gamma',
        with a grid=[1, 1, 1] and shifts=[0, 0, 0]

        :param kmode: The options are 'cartesian', 'reduced', 'gamma', 'monkhorst-pack'
        :param kpoints_list: Explicit list of k-points (Only for 'cartesian' and 'reduced')
        :param weights: List of weights associated to each k-point (Only for 'cartesian' and 'reduced')
        :param grid: Number of kpoints on each direction (Only for 'gamma' and 'monkhorst-pack')
        :param shifts: Shift applied to the grid (Only for 'gamma' and 'monkhorst-pack')
        :param kvertices: List of vertices for a path of k-points (Only for 'path')
        :param intermediates: Number of intermediates between two vertices (Only for 'path')

        :rtype : Kpoints object
        """
        if kmode.lower() not in ['cartesian', 'reduced', 'gamma', 'monkhorst-pack', 'path']:
            raise ValueError("kmode must be 'cartesian', 'reduced', 'gamma', 'monkhorst-pack' or 'path'")

        self.kmode = kmode.lower()

        # All the variables initially None
        self.kpoints_list = None
        self.weights = None
        self.grid = None
        self.shifts = None
        self.kvertices = None
        self.intermediates = None

        # Explicit list of k-points
        if self.kmode in ['cartesian', 'reduced']:
            self.set_kpoints_list(kpoints_list, weights)
        # Grid specification
        elif self.kmode in ['gamma', 'monkhorst-pack']:
            self.set_grid(grid, shifts)
        # Path along the K-space
        elif self.kmode == 'path':
            if kvertices is None:
                self.kvertices = np.array([[0, 0, 0]])
            else:
                if np.array(kvertices).shape == (3,):
                    self.kvertices = np.array(kvertices).reshape((-1, 3))
                elif np.array(kvertices).shape[1] == 3:
                    self.kvertices = np.array(kvertices)
                else:
                    raise ValueError("Wrong value for kvertices, should be an array with shape (nkpt,3)")
            if intermediates is None:
                self.intermediates = 1
            else:
                self.intermediates = int(intermediates)

    def add_kpt(self, kpoint, weight=1):
        """
        Add one extra kpoint to the kpoints_list in the case of
        kmode 'cartesian' or 'reduced'

        :param kpoint: (list, numpy.array) The kpoint to be added
        :param weight: (int, float) the weight associated to that point
        """
        kpoint = list(kpoint)
        assert len(kpoint) == 3
        assert (self.kmode in ['cartesian', 'reduced'])
        if kpoint not in self.kpoints_list:
            self.kpoints_list.append(kpoint)
            self.weights = np.append(self.weights, weight)
        else:
            position = self.kpoints_list.index(kpoint)
            self.weights[position] = weight

    def set_kpoints_list(self, kpoints_list, weights=None):
        """
        Set an explicit list of kpoints with a proper check of correct arguments

        :param kpoints_list: Explicit list of k-points (Only for 'cartesian' and 'reduced')
        :param weights: List of weights associated to each k-point (Only for 'cartesian' and 'reduced')
        """
        assert (self.kmode in ['cartesian', 'reduced'])
        if kpoints_list is None:
            self.kpoints_list = [[0, 0, 0]]
        else:
            if np.array(kpoints_list).shape == (3,):
                self.kpoints_list = np.array(kpoints_list).reshape((-1, 3))
            elif np.array(kpoints_list).shape[1] == 3:
                self.kpoints_list = np.array(kpoints_list)
            else:
                raise ValueError("Wrong value for kpoints_list, should be an array with shape (nkpt,3)")
        self.kpoints_list = generic_serializer(self.kpoints_list)
        nkpt = len(self.kpoints_list)
        if weights is None:
            self.weights = np.ones(nkpt)
        elif len(np.array(weights)) == nkpt:
            self.weights = np.array(weights)
        else:
            raise ValueError("Wrong value for weights, should be an array with shape (nkpt,)")
        self.weights = list(self.weights)

    def set_grid(self, grid, shifts=None):
        """
        Set a grid of kpoints with a proper check of correct arguments

        :param grid: Number of kpoints on each direction (Only for 'gamma' and 'monkhorst-pack')
        :param shifts: Shift applied to the grid (Only for 'gamma' and 'monkhorst-pack')
        """

        assert (self.kmode in ['gamma', 'monkhorst-pack'])
        if grid is None:
            self.grid = [1, 1, 1]
        elif len(grid) == 3:
            self.grid = [int(x) for x in grid]
        else:
            raise ValueError("Wrong value for grid, should be an array with shape (-1,3)")
        if shifts is None:
            self.shifts = [0, 0, 0]
        else:
            self.shifts = generic_serializer(np.array(shifts).reshape(3))

    @property
    def number_of_kpoints(self):
        return self.nkpt

    @property
    def nkpt(self):
        """
        Number of kpoints
        For a grid the value is the product of the number of kpoints
        on each direction.
        """
        if self.kmode in ['cartesian', 'reduced']:
            return len(self.kpoints_list)
        elif self.kmode in ['gamma', 'monkhorst-pack']:
            return np.prod(self.grid)
        elif self.kmode == 'path':
            return len(self.kvertices) * self.intermediates

    def __str__(self):
        """
        String representation of the KPoints object
        """
        kp = ' Mode  : ' + self.kmode + '\n'
        if self.kmode in ['cartesian', 'reduced']:
            kp += ' Number of k-points: ' + str(self.nkpt) + '\n\n'
            for i in range(self.nkpt):
                kp += (" %10.5f %10.5f %10.5f %15.5f\n"
                       % (self.kpoints_list[i][0], self.kpoints_list[i][1], self.kpoints_list[i][2], self.weights[i]))
        elif self.kmode in ['gamma', 'monkhorst-pack']:
            kp += ' Grid  : ' + str(self.grid) + '\n'
            kp += ' Shift : ' + str(self.shifts) + '\n'
        elif self.kmode == 'path':
            kp += ' Number of intermediates: ' + str(self.intermediates) + '\n\n'
            for vertex in self.kvertices:
                kp += " %10.5f %10.5f %10.5f\n" % (vertex[0], vertex[1], vertex[2])
        return kp

    @property
    def to_dict(self):

        ret = {}
        for i in list(self.__dict__.keys()):
            if eval('self.' + i) is not None:
                ret[i] = eval('self.' + i)
        return ret

    @classmethod
    def from_dict(cls, json_dict):
        if json_dict['kmode'] in ['cartesian', 'reduced']:
            return cls(kmode=json_dict['kmode'], kpoints_list=json_dict['kpoints_list'],
                       weights=json_dict['weights'])
        elif json_dict['kmode'] in ['gamma', 'monkhorst-pack']:
            return cls(kmode=json_dict['kmode'], grid=json_dict['grid'], shifts=json_dict['shifts'])
        elif json_dict['kmode'] == 'path':
            return cls(kmode=json_dict['kmode'], kvertices=json_dict['kvertices'],
                       intermediates=json_dict['intermediates'])

    @staticmethod
    def max_kpoints(structure, kpoints_per_atom=1000):
        """
        Returns the maximal grid of kpoints allowed during the
        kpoints convergence test. Based on the routine used in
        pymatgen.

        :param structure: (pychemia.Structure) The structure to compute the maximal grid
        :param kpoints_per_atom: (int) number of kpoints per atom
        :return:
        """

        lengths = [sqrt(sum(map(lambda y: y ** 2, structure.cell[i]))) for i in range(3)]
        ngrid = kpoints_per_atom / structure.natom
        mult = (ngrid * lengths[0] * lengths[1] * lengths[2]) ** (1.0 / 3.0)
        num_div = [int(ceil(1.0 / lengths[i] * mult)) for i in range(3)]
        num_div = [i if i > 0 else 1 for i in num_div]
        return num_div[0], num_div[1], num_div[2]

    def get_density_of_kpoints(self, lattice):

        assert (self.kmode in ['gamma', 'monkhorst-pack'])
        vol = 1.0 / lattice.volume
        return self.number_of_kpoints / vol

    @classmethod
    def optimized_grid(cls, lattice, kp_density=1E5, nkpoints=None, force_odd=False):
        """
        Returns a KPoints object with kmode='gamma' with a grid optimized for a given
         density or number of kpoints.
         If force_odd is True the grid will be such that there are odd numbers on each direction
         such that 'gamma' is always included on the grid.

        :param lattice: (pychemia.Lattice) A PyChemia Lattice object, if the argument is a pychemia.Structure
                        object, its lattice will be used instead.
        :param kp_density: (float) A density of kpoints (default: 1E5)
        :param nkpoints: (int) Aproximate number of kpoints to use, if this argument is present the kp_density
                                will be ignored
        :param force_odd: (bool) If True, force the grid to be odd for the three directions.
        """
        if nkpoints is not None:
            kp_density = None

        # Reciprocal lattice object
        rlattice = lattice.reciprocal()

        if nkpoints is not None:
            vol = 1.0 / lattice.volume
            kp_density = nkpoints / vol

        factor = int(kp_density ** 1.0 / 3.0)

        while True:
            ini_density = kp_density

            grid = np.array([int(max(ceil(factor * rlattice.a), 1)),
                             int(max(ceil(factor * rlattice.b), 1)),
                             int(max(ceil(factor * rlattice.c), 1))])

            # The volume in reciprocal space is the inverse of the volume in real space
            fin_density = float(np.prod(grid) * lattice.volume)

            if fin_density > ini_density and factor > 0:
                factor -= 1
            else:
                break

        if force_odd:
            for i in range(3):
                if grid[i] % 2 == 0:
                    grid[i] += 1
        return cls(kmode='gamma', grid=grid)
