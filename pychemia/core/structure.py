"""
Definition of the class Structure
This class defines methods to create and manipulate
atomic structures such as molecules, clusters and crystals
"""

import numpy as np
import json
import sys
from math import sin, cos
from itertools import combinations
import itertools

from pychemia import log
from pychemia.core.lattice import Lattice
from pychemia.core.delaunay import get_reduced_bases
from pychemia.core.composition import Composition
from pychemia.utils.computing import unicode2string
from pychemia.utils.periodic import mass, atomic_number, covalent_radius, valence, atomic_symbols
from pychemia.utils.mathematics import matrix_from_eig, vector_set_perpendicular, unit_vector
from multiprocessing import Pool


__author__ = "Guillermo Avendano-Franco"
__copyright__ = "Copyright 2014"
__version__ = "0.1"
__maintainer__ = "Guillermo Avendano-Franco"
__email__ = "gtux.gaf@gmail.com"
__status__ = "Development"
__date__ = "June 10, 2014"


class Structure():
    """
    Define an object that contains information about atomic positions,
    cell parameters and periodicity and provides methods to manipulate
    those elements

    Represents a molecule, cluster, wire, slab or crystal structure
    The positions of the atoms and their atomic symbols are declared
    in 'positions' and 'symbols' respectively.

    For periodic structures, the 'periodicity' can be declared.
    and cell parameters in 'cell'

    Magnetic moments can be associated in the array vector_info['magnetic_moments'].

    This object contains no dynamical information. That information
    is supported by the child class DynamicStructure

    """

    def __init__(self, **kwargs):
        """
        :param comment:
        :param natom:
        :param symbols:
        :param periodicity:
        :param cell:
        :param positions:
        :param reduced:
        :param mag_moments:
        :param kwargs:

        Args:

        natom      : Integer with number of atoms
        symbols    : String list of atom symbols
        positions  : Array of atomic positions
                     Those are dimensional units
        cell       : Dimensional cell (3x3 matrix)
        reduced    : Cell-scaled positions
                     Dimensionless
        mag_moments: Array of Magnetic moments
        periodicity: Periodicity on each direction
                     3-array of booleans
        symm   : Symmetry information like point symm operations
                     and space group
        name       : Free text to identify structure (Only one line, max 50 chars)
        comment    : Free text to identify structure

        Examples:

>>> import pychemia
>>> a = pychemia.Structure()
>>> print a
Empty structure
>>> a = pychemia.Structure(symbols=['Xe'])
>>> print a.natom
1
>>> d = 1.104
>>> a = pychemia.Structure(symbols=['N', 'N'], positions=[[0, 0, 0], [0, 0, d]], periodicity=False)
>>> print a.natom
2
>>> a = 4.05
>>> b = a/2
>>> fcc = pychemia.Structure(symbols=['Au'], cell=[[0, b, b], [b, 0, b], [b, b, 0]], periodicity=True)
>>> print fcc.natom
1
        """
        self.vector_info = {}
        self.name = None
        self.comment = None
        self.natom = None
        self.symbols = None
        self.positions = None
        self.reduced = None
        self.cell = None
        self.periodicity = None
        self.vector_info['mag_moments'] = None
        self.sites = None
        self.occupancies = None

        self._lattice = None
        self._composition = None

        # Fill the values from args
        if 'name' in kwargs and kwargs['name'] is not None:
            self.name = kwargs['name'].split('\n')[0][:50]
        if 'comment' in kwargs:
            self.comment = kwargs['comment']
        if 'natom' in kwargs:
            self.natom = int(kwargs['natom'])
        if 'symbols' in kwargs:
            self.symbols = kwargs['symbols']
        if 'periodicity' in kwargs:
            periodicity = kwargs['periodicity']
            self.set_periodicity(periodicity)
        if 'cell' in kwargs:
            cell = np.array(kwargs['cell'])
            self.set_cell(cell)
        if 'positions' in kwargs:
            positions = np.array(kwargs['positions'])
            self.set_positions(positions)
        if 'reduced' in kwargs:
            reduced = np.array(kwargs['reduced'])
            self.set_reduced(reduced)
        if 'mag_moments' in kwargs:
            self.set_mag_moments(np.array(kwargs['mag_moments']))
        if 'occupancies' in kwargs:
            self.occupancies = kwargs['occupancies']
        if 'sites' in kwargs:
            self.sites = kwargs['sites']

        # Lets autocomplete the missing information
        self._autocomplete()

        if not self._check():
            print('Arguments non consistent')

    def __len__(self):
        return self.natom

    def __str__(self):
        if self.natom == 0:
            xyz = 'Empty structure'
        else:
            xyz = str(self.natom) + '\n\n'
            if self.is_crystal:
                xyz += ' Symb  (             Positions            ) [     Cell-reduced coordinates     ]\n'
            else:
                xyz += ' Symb  (             Positions            )\n'

            for i in range(self.natom):
                if self.is_crystal:
                    xyz += (" %4s  ( %10.4f %10.4f %10.4f ) [ %10.4f %10.4f %10.4f ]\n"
                            % (self.symbols[i],
                               self.positions[i, 0],
                               self.positions[i, 1],
                               self.positions[i, 2],
                               self.reduced[i, 0],
                               self.reduced[i, 1],
                               self.reduced[i, 2]))
                else:
                    xyz += (" %4s  ( %10.4f %10.4f %10.4f )\n"
                            % (self.symbols[i],
                               self.positions[i, 0],
                               self.positions[i, 1],
                               self.positions[i, 2]))

            if self.periodicity[0] or self.periodicity[1] or self.periodicity[2]:
                xyz += '\nPeriodicity: '
                if self.periodicity[0]:
                    xyz += ' X'
                if self.periodicity[1]:
                    xyz += ' Y'
                if self.periodicity[2]:
                    xyz += ' Z'
                xyz += '\n\nLattice vectors:\n'
                for i in range(3):
                    xyz += (" %10.4f %10.4f %10.4f\n"
                            % (self.cell[i, 0], self.cell[i, 1], self.cell[i, 2]))
            else:
                xyz += '\nNon-periodic structure'

        return xyz

    def __repr__(self):
        ret = 'Structure(symbols=' + str(self.symbols)
        if self.is_periodic:
            if np.all(np.diag(self.cell.diagonal()) == self.cell):
                if np.max(self.cell.diagonal()) == np.min(self.cell.diagonal()):
                    ret += ', cell=' + str(self.cell[0, 0])
                else:
                    ret += ', cell=' + str(self.cell.diagonal().tolist())
            else:
                ret += ', cell=' + str(self.cell.tolist())
        ret += ', positions=' + str(self.positions.tolist())
        if all([self.periodicity[0] == item for item in self.periodicity]):
            ret += ', periodicity=' + str(self.periodicity[0])
        else:
            ret += ', periodicity=' + str(self.periodicity)
        ret += ')'
        return ret

    def __iter__(self):
        return iter(SiteSet(self))

    def _autocomplete(self):
        if self.natom is None:
            if self.positions is not None:
                self.natom = len(self.positions)
            elif self.reduced is not None:
                self.natom = len(self.reduced)
            elif self.symbols is not None:
                self.natom = len(self.symbols)
            else:
                self.natom = 0

        if self.symbols is None and self.natom == 0:
            self.symbols = []

        if self.periodicity is None:
            self.set_periodicity(True)

        if self.cell is None and self.is_periodic:
            self.set_cell(1)

        if self.positions is None:
            if self.reduced is not None:
                self.reduced2positions()
            else:
                if self.natom == 0:
                    self.positions = np.array([])
                elif self.natom == 1:
                    self.positions = np.array([[0.0, 0.0, 0.0]])
                else:
                    raise ValueError('Positions must be present for more than 1 atom')

        if self.reduced is None and self.is_crystal:
            if self.positions is not None and self.natom > 0:
                self.positions2reduced()
            else:
                self.reduced = np.array([])

        if self.sites is None:
            self.sites = range(self.natom)

        if self.occupancies is None:
            self.occupancies = np.ones(self.natom)

    def _check(self):
        check = True

        if len(self.symbols) != self.natom:
            print('Error: Bad symbols')
            check = False
        if len(self.positions) != self.natom:
            print('Error: Bad positions')
            check = False
        if self.is_crystal and len(self.reduced) != self.natom:
            print('Error: Bad reduced')
            check = False
        if self.vector_info['mag_moments'] is not None and len(self.vector_info['mag_moments']) != self.natom:
            print('Error: Bad mag_moments')
            check = False

        return check

    def add_atom(self, name, coordinates, option='cartesian'):
        """
        Add an atom with a given 'name' and cartesian or reduced 'position'
        The atom will be added at the end of the list of atoms in the Structure

        :param name: (str)
        :param coordinates: (list, numpy.array)
        :param option: (str)
        """
        assert (name in atomic_symbols)
        assert (option in ['cartesian', 'reduced'])
        self.symbols.append(name)
        self.natom += 1
        self._composition = None

        if option == 'cartesian':
            if self.natom == 0:
                self.positions = np.array(coordinates).reshape([-1, 3])
            else:
                self.positions = np.append(self.positions, coordinates).reshape([-1, 3])
            self.positions2reduced()
        elif option == 'reduced':
            if self.natom == 0:
                self.reduced = np.array(coordinates).reshape([-1, 3])
            else:
                self.reduced = np.append(self.reduced, coordinates).reshape([-1, 3])
            self.reduced2positions()

    def del_atom(self, index):
        """
        Removes the atom with the given index

        :param index:

        :return:
        """
        assert (abs(index) < self.natom)
        self.symbols.pop(index)
        np.delete(self.positions, index, 0)
        np.delete(self.reduced, index, 0)
        self.natom -= 1
        self._composition = None

    def center_mass(self, list_of_atoms=None):
        """
        Computes the center of mass (CM) of the XYZ object or
        a partial list of atoms. The default is to compute the
        CM of all the atoms in the object, if a list
        is enter only those in the list will be included for the CM
        Return the CM as a numpy array
        """
        if list_of_atoms is None:
            list_of_atoms = range(self.natom)

        total_mass = 0.0
        center_of_mass = np.zeros(3)
        if self.natom == 0:
            return center_of_mass

        atomicnumber = atomic_number(list(self.symbols))

        for i in range(self.natom):
            if i in list_of_atoms:
                total_mass = total_mass + mass(atomicnumber[i])
                center_of_mass = center_of_mass + mass(atomicnumber[i]) * self.positions[i]

        return center_of_mass / total_mass

    def rotation(self, tx, ty, tz):
        """
        Rotate the molecule in the three directions
        """

        rotationx = np.array([[1, 0, 0], [0, cos(tx), -sin(tx)], [0, sin(tx), cos(tx)]])
        rotationy = np.array([[cos(ty), 0, sin(ty)], [0, 1, 0], [-sin(ty), 0, cos(ty)]])
        rotationz = np.array([[cos(tz), -sin(tz), 0], [sin(tz), cos(tz), 0], [0, 0, 1]])

        rotation = np.dot(np.dot(rotationx, rotationy), rotationz)

        for i in range(self.natom):
            self.positions[i] = np.dot(rotation, self.positions[i])

    def get_cell(self):
        if self._lattice is None:
            self._lattice = Lattice(self.cell)
        return self._lattice

    @property
    def lattice(self):
        return self.get_cell()

    def get_composition(self, gcd=True):
        """
        Computes the composition of the Structure
        as the count of each species in the cell
        If gcd is True the values are divided by the 
        greatest common divisor 

        :param gcd: bool

        :rtype : Composition
        """
        if self._composition is None:
            species = {}
            for atom in self.symbols:
                if atom in species:
                    species[atom] += 1
                else:
                    species[atom] = 1
            self._composition = Composition(species)
        return self._composition

    def positions2reduced(self):
        """
        Computes the cell-reduced coordinates from the
        cartesian dimensional coordinates
        """
        self.reduced = np.linalg.solve(self.cell.T, self.positions.T).T
        for i in range(3):
            if self.periodicity[i]:
                self.reduced[:, i] %= 1.0

    def reduced2positions(self):
        """
        Computes the dimensional cartesian coordinates
        from the adimensional cell-reduced coordinates
        """
        self.positions = np.dot(self.reduced, self.cell)

    def relocate_to_cm(self, list_of_atoms=None):
        """
        Relocates the system of atoms to the center of mass
        a partial list of atoms can be used to compute
        the center, but all the atoms are moved to the
        computed center
        """
        cm = self.center_mass(list_of_atoms)
        self.positions = self.positions - cm

    def get_distance(self, iatom, jatom, with_periodicity=True, tolerance=1e-5):
        """
        Calculates the distance between 2 atom, identified by index
        iatom and jatom

        :param iatom: (int) index of first atom
        :param jatom: (int) index of second atom
        :param with_periodicity: (bool) if the periodic images should be considered to compute the shortest distance
        :param tolerance: (float) Tolerance for the bases reduction

        :rtype : (float) distance between iatom and jatom
        """

        if with_periodicity:
            reduced_bases = get_reduced_bases(self.cell, tolerance)
            scaled_pos = np.dot(self.positions, np.linalg.inv(reduced_bases))
            # move scaled atomic positions into -0.5 < r <= 0.5
            for pos in scaled_pos:
                pos -= pos.round()

            # Look for the shortest one in surrounded 3x3x3 cells
            distances_list = []
            for i in (-1, 0, 1):
                for j in (-1, 0, 1):
                    for k in (-1, 0, 1):
                        distances_list.append(np.linalg.norm(
                            np.dot(scaled_pos[iatom] - scaled_pos[jatom] +
                                   np.array([i, j, k]), reduced_bases)))
            ret = min(distances_list)

        else:
            posi = self.positions[iatom]
            posj = self.positions[jatom]
            ret = np.linalg.norm(posi - posj)

        return ret

    @staticmethod
    def random_cell(composition, method='stretching', stabilization_number=20, nparal=5):
        """
        Generate a random cell
        There are two algorithms implemented:

        scaling: Generate a random cell and random distribution of atoms and
                scale the lattice to separate the atoms.

        stretching: Generating a random cell and random distribution of atoms
                    and stretching their bonds until the distance between any
                    two atoms is always greater than the sum of covalent radius.
        """
        comp = Composition(composition)
        log.debug('Generating a random structure with composition: ' + str(comp.composition))
        natom = comp.natom
        symbols = comp.symbols

        best_volume = sys.float_info.max
        best_volume = float('inf')
        best_structure = None
        optimal_volume = comp.covalent_volume('cubes')
        stabilization_history = 0
        pool = Pool(processes=nparal)

        trial = 0
        while stabilization_history < stabilization_number:
            args = list(best_volume * np.ones(10))

            ret = pool.map(worker_star, itertools.izip(itertools.repeat(method), itertools.repeat(composition), args))

            ngood = 0
            for structure in ret:
                if structure is not None:
                    ngood += 1
                    if best_structure is not None:
                        if structure.volume < best_structure.volume:
                            best_structure = structure
                    else:
                        best_structure = structure

            # log.debug('Good structures: %d/10  Best volume: %7.3f' % (ngood, best_structure.volume))
            if best_structure is not None and best_volume > best_structure.volume:
                best_volume = best_structure.volume
                stabilization_history = 0
            else:
                stabilization_history += 1
            trial += 1

            # log.debug('Trial: %4d  Volume: %7.2f  Optimal Volume: %7.2f  Ratio: %5.2f' %
            #          (trial, best_volume, optimal_volume, best_volume/optimal_volume))

        pool.close()

        if best_structure is not None:
            # Analysis of the quality for the best structure
            rpos = best_structure.reduced
            for i, j in combinations(range(natom), 2):
                distance = best_structure.lattice.minimal_distance(rpos[i], rpos[j])
                covalent_distance = sum(covalent_radius([symbols[i], symbols[j]]))
                if distance < covalent_distance:
                    log.debug('Covalent distance: %7.4f  Minimal distance: %7.4f  Difference: %7.3e' %
                              (covalent_distance, distance, covalent_distance - distance))

        return best_structure

    def adjust_reduced(self):
        for i in range(self.natom):
            for j in range(3):
                for value in [0.5, 0.25, 0.75, 0.125]:
                    if abs(value - self.reduced[i, j]) < 1E-4:
                        self.reduced[i, j] = value
        self.reduced2positions()

    def set_cell(self, cell):
        """
        Set the vectors defining the cell

        Args:
            cell: A matrix with the 3 unit cell
            vectors
        """
        npcell = np.array(cell)
        if npcell.shape == () or npcell.shape == (1,):
            self.cell = npcell * np.eye(3)
        elif npcell.shape == (3,):
            self.cell = np.diag(npcell)
        else:
            self.cell = np.array(cell).reshape([3, 3])
        self._lattice = None

    def set_mag_moments(self, mag_moments):
        """
        Set the magnetic moments with one vector on each
        atom

        Args:
            mag_moments: List or numpy array with one
            vector for each atom.
            The values will be converted into a numpy array
        """
        self.vector_info['mag_moments'] = np.array(mag_moments).reshape([-1, 3])

    def set_periodicity(self, periodicity):
        """
        Set periodicity of the structure

        Args:
            periodicity: (Boolean) a single value
            means that the structure has that
            periodicity all along the 3 directions.
            Otherwise a list of 3 booleans is
            required
        """
        if isinstance(periodicity, bool):
            self.periodicity = 3 * [periodicity]
        elif isinstance(periodicity, list) and len(periodicity) == 1:
            self.periodicity = 3 * periodicity
        else:
            self.periodicity = list(periodicity)

    def set_positions(self, positions):
        """
        Set the positions of the atoms
        This contains dimensional values
        in cartesian coordinates

        Args:
            positions: A array of 3 vectors
            with dimensional coordinates
        """
        self.positions = np.array(positions).reshape([-1, 3])

    def set_reduced(self, reduced):
        """
        Set the reduced positions of the atoms
        This contains adimensional values
        relative to cell vectors

        Args:
            positions: A array of 3 vectors
            with adimensional coordinates
        """
        self.reduced = np.array(reduced).reshape([-1, 3])

    def sort_sites_using_list(self, sorted_indices):
        self.symbols = list(np.array(self.symbols)[sorted_indices])
        self.positions = self.positions[sorted_indices]
        self.reduced = self.reduced[sorted_indices]
        if self.vector_info is not None:
            for vi in self.vector_info:
                if self.vector_info[vi] is not None:
                    self.vector_info[vi] = self.vector_info[vi][sorted_indices]

    def sort_sites(self):
        # Only implemented for perfect structures
        assert self.is_perfect

        # First: Sort sites using the distance to the origin
        sorted_indices = np.array([np.linalg.norm(self.positions[i]) for i in range(self.nsites)]).argsort()
        self.sort_sites_using_list(sorted_indices)

        # Second: Sort again using the atomic number
        sorted_indices = np.array([atomic_number(x) for x in self.symbols]).argsort()
        self.sort_sites_using_list(sorted_indices)

        # The sites are ordered by atomic number first and distance to the origin second

    def sort_axes(self):
        """
        Sort the lattice vectors in decremental order of their size.
        'a' will be the longest lattice vector
        'c' he shortest
        """

        sorted_indices = self.lattice.lengths.argsort()[::-1]
        self.set_cell(self.cell[sorted_indices])
        self.reduced = self.reduced[:, sorted_indices]

    def align_with_axis(self, axis=0, round_decimals=8):
        lattice = self.lattice
        lattice.align_with_axis(axis=axis, round_decimals=round_decimals)
        self.set_cell(lattice.cell)
        self.reduced2positions()

    def canonical_form(self):

        self.sort_sites()
        self.sort_axes()
        self.align_with_axis()

    def supercell(self, size):
        """
        Creates a supercell, replicating the positions
        of atoms in the x,y,z directions a number of
        size=(nx,ny,nz) times
        """
        new_natom = np.prod(size) * self.natom
        new_symbols = np.prod(size) * self.symbols
        new_positions = np.zeros((new_natom, 3))

        index = 0
        for i in range(size[0]):
            for j in range(size[1]):
                for k in range(size[2]):
                    for n in range(self.natom):
                        new_positions[index] = self.positions[n] + (
                            i * self.cell[0] + j * self.cell[1] + k * self.cell[2])
                        index += 1
        new_cell = np.zeros((3, 3))
        new_cell[0] = size[0] * self.cell[0]
        new_cell[1] = size[1] * self.cell[1]
        new_cell[2] = size[2] * self.cell[2]
        return Structure(symbols=new_symbols, positions=new_positions, cell=new_cell)

    def copy(self):
        """
        Get a copy of the object
        """
        copy_struct = Structure(name=self.name, comment=self.comment, natom=self.natom, symbols=self.symbols,
                                periodicity=self.periodicity, cell=self.cell, positions=self.positions,
                                reduced=self.reduced, vector_info=self.vector_info, sites=self.sites,
                                occupancies=self.occupancies)
        return copy_struct

    def plot(self, figname='None', size=(300, 325), save=False):
        from mayavi import mlab

        fig = mlab.figure(size=size)
        figure = mlab.gcf()
        fig.scene.disable_render = True
        figure.scene.background = (0.0, 0.0, 0.0)
        mlab.view(30, 30)
        assert (self.natom > 0)

        x = self.positions[:, 0]
        y = self.positions[:, 1]
        z = self.positions[:, 2]
        cr = covalent_radius(self.symbols)

        mlab.points3d(x, y, z, cr, scale_factor=1)

        if self.is_crystal:
            frame, line1, line2, line3 = self.get_cell().get_path()

            mlab.plot3d(frame[:, 0], frame[:, 1], frame[:, 2], tube_radius=.02, color=(1, 1, 1))
            mlab.plot3d(line1[:, 0], line1[:, 1], line1[:, 2], tube_radius=.02, color=(1, 1, 1))
            mlab.plot3d(line2[:, 0], line2[:, 1], line2[:, 2], tube_radius=.02, color=(1, 1, 1))
            mlab.plot3d(line3[:, 0], line3[:, 1], line3[:, 2], tube_radius=.02, color=(1, 1, 1))

        fig.scene.disable_render = False
        if save:
            mlab.savefig(figname)
        return figure

    def to_dict(self):

        ret = {'name': self.name,
               'comment': self.comment,
               'natom': self.natom,
               'symbols': self.symbols,
               'periodicity': self.periodicity,
               'cell': self.cell.tolist(),
               'positions': self.positions.tolist(),
               'reduced': self.reduced.tolist(),
               'vector_info': self.vector_info,
               'nspecies': len(self.species),
               'density': self.density,
               'formula': self.formula,
               'sites': list(self.sites),
               'occupancies': list(self.occupancies)}
        return ret

    @staticmethod
    def from_dict(structdict):

        name = structdict['name']
        comment = structdict['comment']
        natom = structdict['natom']
        symbols = unicode2string(structdict['symbols'])
        periodicity = structdict['periodicity']
        cell = np.array(structdict['cell'])
        positions = np.array(structdict['positions'])
        reduced = np.array(structdict['reduced'])
        vector_info = structdict['vector_info']
        if 'sites' in structdict:
            sites = structdict['sites']
        else:
            sites = range(natom)
        if 'occupancies' in structdict:
            occupancies = structdict['occupancies']
        else:
            occupancies = list(np.ones(natom))
        return Structure(name=name, comment=comment, natom=natom, symbols=symbols, periodicity=periodicity, cell=cell,
                         positions=positions, reduced=reduced, vector_info=vector_info, sites=sites,
                         occupancies=occupancies)

    def save_json(self, filename):

        filep = open(filename, 'w')
        json.dump(self.to_dict(), filep, sort_keys=True, indent=4, separators=(',', ': '))
        filep.close()

    @staticmethod
    def load_json(filename):

        filep = open(filename, 'r')
        structdict = unicode2string(json.load(filep))
        filep.close()
        return Structure.from_dict(structdict)

    def distance2(self, atom1, atom2):
        assert (isinstance(atom1, int))
        assert (isinstance(atom2, int))
        assert (atom1 < self.natom)
        assert (atom2 < self.natom)
        return self.lattice.distance2(self.reduced[atom1], self.reduced[atom2])

    def valence_electrons(self):
        ret = 0
        for key, value in self.composition.items():
            ret += value * valence(key)
        return ret

    def __eq__(self, other):
        if self.natom != other.natom:
            ret = False
        elif not np.array_equal(self.cell, other.cell):
            ret = False
        elif not np.array_equal(self.positions, other.positions):
            ret = False
        elif not np.array_equal(self.reduced, other.reduced):
            ret = False
        elif not np.array_equal(self.periodicity, other.periodicity):
            ret = False
        else:
            ret = True
        return ret

    def __ne__(self, other):
        return not self.__eq__(other)

    @property
    def is_perfect(self):
        """
        Return True if two conditions are met:

        1. The number of sites is equal to
        the number of atoms. ie there is no more than
        one atom on each site.

        2. All the occupancies are equal to 1

        :rtype : bool

        :return: bool
        """
        return self.natom == self.nsites and min(self.occupancies) == 1

    @property
    def is_periodic(self):
        """
        Return True if the Structure is periodic in any direction
        False for non-periodic structures

        :rtype : bool

        :return: bool
        """
        return any(self.periodicity)

    @property
    def is_crystal(self):
        """
        True if structure is periodic in all directions
        False otherwise

        :rtype : bool

        :return: bool
        """
        if not self.is_periodic:
            return False
        else:
            return self.get_cell().periodic_dimensions == 3

    @property
    def composition(self):
        """
        Dictionary with the composition, the keys are the species and the values
        represent the number of atoms of that specie

        :rtype : dict

        :return: dict
        """
        return self.get_composition().composition

    @property
    def formula(self):
        """
        String with the chemical formula

        :rtype: str

        :return: str
        """
        return self.get_composition().formula

    @property
    def density(self):
        """
        Computes the density of the cell

        :rtype: float

        :return: float
        """
        return sum(np.array(mass(self.symbols))) / self.volume

    @property
    def volume(self):
        """
        Computes the volume of the cell

        :rtype: float

        :return: float
        """
        return abs(np.linalg.det(self.cell))

    @property
    def species(self):
        return self.get_composition().species

    @property
    def nspecies(self):
        return len(self.get_composition().species)

    @property
    def nsites(self):
        return len(self.positions)


def load_structure_json(filename):
    ret = Structure()
    ret.load_json(filename)
    return ret


class SiteSet():
    def __init__(self, structure):

        self.structure = structure
        self.sitelist = []

        for isite in range(structure.nsites):
            if structure.sites.count(isite) > 1:
                symbols = []
                occupancies = []
                for jatom in range(structure.natom):
                    if structure.sites[jatom] == isite:
                        symbols.append(structure.symbols[jatom])
                        occupancies.append(structure.occupancies[jatom])
                position = structure.positions[isite]
                reduced = structure.reduced[isite]
            else:
                symbols = [structure.symbols[isite]]
                occupancies = [structure.occupancies[isite]]
                position = structure.positions[isite]
                reduced = structure.reduced[isite]
            self.sitelist.append(Site(symbols=symbols, occupancies=occupancies, position=position, reduced=reduced))

    def __iter__(self):
        return iter(self.sitelist)


class Site():
    def __init__(self, symbols, occupancies, position, reduced):

        if isinstance(symbols, list):
            self.symbols = symbols
        else:
            self.symbols = [symbols]

        if isinstance(occupancies, list):
            self.occupancies = occupancies
        else:
            self.occupancies = [occupancies]

        assert (len(self.occupancies) == len(self.symbols))

        self.position = position
        self.reduced = reduced

    def __repr__(self):
        ret = 'Site(symbols=' + repr(self.symbols)
        ret += ',occupancies=' + repr(self.occupancies)
        ret += ',position=' + repr(self.position)
        ret += ',reduced=' + repr(self.reduced)
        ret += ')'
        return ret

    def __str__(self):
        return repr(self)


def worker_star(x):
    return random_structure(*x)


def random_structure(method, composition, best_volume=1E10):

    comp = Composition(composition)
    natom = comp.natom
    symbols = comp.symbols

    new_structure = None
    if method == 'scaling':
        lattice = Lattice.random_cell(comp)
        # Random reduced positions
        rpos = np.random.rand(natom, 3)
        mins = [min(rpos[:, i]) for i in range(3)]
        rpos -= mins

        factor = 1.0
        for i in range(len(rpos)):
            covalent_dim = 2.0 * covalent_radius(symbols[i])

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
                covalent_dim = sum(covalent_radius([symbols[i], symbols[j]]))
                this_factor = covalent_dim / distance
                if this_factor > factor:
                    factor = this_factor
        a = lattice.a
        b = lattice.b
        c = lattice.c
        alpha = lattice.alpha
        beta = lattice.beta
        gamma = lattice.gamma
        lattice = Lattice().from_parameters_to_cell(factor * a, factor * b, factor * c, alpha, beta, gamma)

    elif method == 'stretching':

        lattice = Lattice.random_cell(comp)
        # Random reduced positions
        rpos = np.random.rand(natom, 3)
        mins = [min(rpos[:, i]) for i in range(3)]
        rpos -= mins
        # Dummy array that always will be overwritten
        eigv = np.array([1, 0, 0])

        for i, j in combinations(range(natom), 2):
            ret = lattice.distance2(rpos[i], rpos[j])
            mindist = sys.float_info.max
            for k in ret:
                if 0 < ret[k]['distance'] < mindist:
                    mindist = ret[k]['distance']
                    eigv = ret[k]['image']

            covalent_distance = sum(covalent_radius([symbols[i], symbols[j]]))
            if mindist < covalent_distance:
                factor = 1.1 * covalent_distance / mindist
                v1, v2, v3 = vector_set_perpendicular(eigv)
                matrixA = matrix_from_eig(v1, v2, v3, factor, 1, 1)
                lattice = Lattice(np.dot(matrixA, lattice.cell))

    if lattice.volume < best_volume:
        test = True
        for i in range(natom):
            for j in range(i + 1, natom):
                distance = lattice.minimal_distance(rpos[i], rpos[j])
                covalent_dim = sum(covalent_radius([symbols[i], symbols[j]]))
                if distance < covalent_dim:
                    test = False
        if test:
            new_structure = Structure(symbols=symbols, reduced=rpos, cell=lattice.cell, periodicity=True)
    return new_structure
    # End of Worker
