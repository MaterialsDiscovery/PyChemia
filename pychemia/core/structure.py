"""
Structure contains the description of a set of atoms in space, for periodic structures it adds a lattice.
"""
try:
    import itertools.izip as zip
except ImportError:
    pass

import json
import os
import struct
import sys
import numpy as np
from collections.abc import MutableSequence
from itertools import combinations, repeat
from math import sin, cos
from multiprocessing import Pool
from pychemia import pcm_log
from pychemia.crystal.lattice import Lattice
from pychemia.core.composition import Composition
from pychemia.core.delaunay import get_reduced_bases
from pychemia.utils.computing import deep_unicode
from pychemia.utils.periodic import mass, atomic_number, covalent_radius, valence, atomic_symbols
import scipy.spatial


class Structure(MutableSequence):
    """
    Define an object that contains information about atomic positions,
    cell parameters and periodicity and provides methods to manipulate
    those elements

    A Structure is basically a set of sites with eventually a lattice
    each site could have one or more species with occupations equal
    or lower than one.

    A Structure could represents a molecule, cluster, wire, slab,
    crystal structure or alloy with define sites.
    The positions of the atoms and their atomic symbols are declared
    in 'positions' and 'symbols' respectively.

    For periodic structures, the 'periodicity' can be declared.
    and cell parameters in 'cell'

    Magnetic moments can be associated in the array vector_info['magnetic_moments'].
    """

    def __init__(self, natom=None, symbols=None, periodicity=False, cell=None, positions=None, reduced=None,
                 mag_moments=None, occupancies=None, sites=None, name=None, comment=None, vector_info=None):
        """ Structure is a container for geometric structure and composition for both periodic and non periodic
        atomic structures

        :param natom: Number of atoms
        :param symbols: List of atomic symbols
        :param periodicity: (True) if the structure is periodic, False for finite structures and a list of booleans for
        structures that are periodic only along specific directions.
        :param cell: The cell parameters, could be scalar for a cubic cell, a list of 3 numbers for a orthogonal cell
        or a complete numpy array or list of lists. Each row is considered a cell vector.
        :param positions: Array of rows with 3 elements for the positions of atoms in cartesian coordinates.
        :param reduced: Positions of atoms as reduced coordinates relative to the cell vectors ie cell-scaled in range
        [0,1]
        :param mag_moments: Magnetic moments for atoms in the structure
        :param occupancies: Atomic occupancies. 1 by default, lower values for vacancies and non-perfect crystals.
        :param sites: Atomic sites

        >>> a = Structure()
        >>> print(a)
        Empty structure
        >>> a = Structure(symbols=['Xe'])
        >>> print(a.natom)
        1
        >>> d = 1.104
        >>> a = Structure(symbols=['N', 'N'], positions=[[0, 0, 0], [0, 0, d]], periodicity=False)
        >>> print(a.natom)
        2
        >>> a = 4.05
        >>> b = a/2
        >>> fcc = Structure(symbols=['Au'], cell=[[0, b, b], [b, 0, b], [b, b, 0]], periodicity=True)
        >>> print(fcc.natom)
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
        self.vector_info = None

        # By default the number of atoms will be the value given or zero except if other information overrules
        # that value
        if natom is not None:
            self.natom = int(natom)
        else:
            self.natom = 0

        if symbols is not None:
            if symbols in atomic_symbols:
                self.symbols = [symbols]
            else:
                for iatom in list(symbols):
                    assert(iatom in atomic_symbols)
                self.symbols = list(symbols)

            self.natom = len(self.symbols)
        else:
            if self.natom != 0:
                raise ValueError('List of atomic symbols not provided for structure with %d atoms', self.natom)

        # No periodicity will be assumed except if cell or reduced coordinates are provided
        if periodicity is None:
            periodicity = 3*[False]
        if isinstance(periodicity, bool):
            periodicity = 3*[periodicity]
        self.set_periodicity(periodicity)

        if cell is not None:
            cell = np.array(cell)
            self.set_cell(cell)
            self.periodicity = 3*[True]

        if positions is not None:
            positions = np.array(positions)
            self.set_positions(positions)
        if reduced is not None:
            reduced = np.array(reduced)
            self.set_reduced(reduced)
            self.periodicity = 3*[True]
        if mag_moments is not None:
            self.set_mag_moments(np.array(mag_moments))
        if occupancies is not None:
            self.occupancies = list(occupancies)
        if sites is not None:
            self.sites = sites

        if vector_info is None:
            self.vector_info = {'mag_moments': None}
        else:
            self.vector_info = vector_info

        self.name = name
        self.comment = comment

        # This routine completes the missing values and makes all the values coherent.
        self._autocomplete()

        if not self._check():
            raise ValueError('Arguments non consistent')

    def __len__(self):
        """ Number of sites in structure.
        In for perfect crystals it will match the number of atoms as each atomic site will have a single atom.

        :return: Number of sites in structure
        :rtype; int

        >>> st = Structure(symbols=['H', 'O'], positions= [[0,0,0], [0,0,1]])
        >>> len(st)
        2
        """
        return self.nsites

    def __str__(self):
        """String representation of Structure

        :return: Human readable text for Structure

        >>> st = Structure(symbols=['Xe'])
        >>> print(st)
        1
        <BLANKLINE>
        Symb  (             Positions            )
          Xe  (     0.0000     0.0000     0.0000 )
        <BLANKLINE>
        Non-periodic structure

        >>> a = 4.05
        >>> b = a/2
        >>> fcc = Structure(symbols=['Au'], cell=[[0, b, b], [b, 0, b], [b, b, 0]], periodicity=True)
        >>> print(fcc)
        1
        <BLANKLINE>
        Symb  (             Positions            ) [     Cell-reduced coordinates     ]
          Au  (     0.0000     0.0000     0.0000 ) [     0.0000     0.0000     0.0000 ]
        <BLANKLINE>
        Periodicity:  X Y Z
        <BLANKLINE>
        Lattice vectors:
             0.0000     2.0250     2.0250
             2.0250     0.0000     2.0250
             2.0250     2.0250     0.0000
        <BLANKLINE>
        """
        if self.natom == 0:
            xyz = 'Empty structure'
        else:
            xyz = str(self.natom) + '\n\n'
            if self.is_crystal:
                xyz += 'Symb  (             Positions            ) [     Cell-reduced coordinates     ]\n'
            else:
                xyz += 'Symb  (             Positions            )\n'

            for i in range(self.natom):
                if self.is_crystal:
                    xyz += ("%4s  ( %10.4f %10.4f %10.4f ) [ %10.4f %10.4f %10.4f ]\n"
                            % (self.symbols[i],
                               self.positions[i, 0],
                               self.positions[i, 1],
                               self.positions[i, 2],
                               self.reduced[i, 0],
                               self.reduced[i, 1],
                               self.reduced[i, 2]))
                else:
                    xyz += ("%4s  ( %10.4f %10.4f %10.4f )\n"
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
        """
        Evaluatable representation of Structure

        :return: String representation of the structure
        :rtype: str

        >>> st1 = Structure(symbols=['H'])
        >>> st2 = eval(repr(st1))
        >>> st1 == st2
        True
        >>> st = Structure(symbols='He', cell=[2,2,2])
        >>> st
        Structure(symbols=['He'], cell=2, reduced=[[0.0, 0.0, 0.0]], periodicity=True)
        """
        ret = 'Structure(symbols=' + str(self.symbols)
        if self.is_periodic:
            if np.all(np.diag(self.cell.diagonal()) == self.cell):
                if np.max(self.cell.diagonal()) == np.min(self.cell.diagonal()):
                    ret += ', cell=' + str(self.cell[0, 0])
                else:
                    ret += ', cell=' + str(self.cell.diagonal().tolist())
            else:
                ret += ', cell=' + str(self.cell.tolist())
            ret += ', reduced=' + str(self.reduced.tolist())
        else:
            ret += ', positions=' + str(self.positions.tolist())
        if all([self.periodicity[0] == item for item in self.periodicity]):
            ret += ', periodicity=' + str(self.periodicity[0])
        else:
            ret += ', periodicity=' + str(self.periodicity)
        ret += ')'
        return ret

    def __delitem__(self, key):
        self.del_atom(key)

    def __setitem__(self, key, value):
        self.add_atom(value['symbols'], value['positions'])

    def __getitem__(self, item):
        return SiteSet(self)[item]

    def __iter__(self):
        return iter(SiteSet(self))

    def insert(self, index, value):
        self.add_atom(value['symbols'], value['positions'])

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
            self.occupancies = self.natom * [1.0]

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
        if self.is_periodic:
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
                total_mass += mass(atomicnumber[i])
                center_of_mass += mass(atomicnumber[i]) * self.positions[i]

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

        :param list_of_atoms: (list) List of atoms that will be considered for computing the center of mass
                              by default all atoms are included
        """
        cm = self.center_mass(list_of_atoms)
        self.positions += -cm

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
    def random_cell(composition, method='stretching', stabilization_number=20, nparal=5, periodic=True,
                    factor_optimal_volume=8):
        """
        Generate a random cell
        There are two algorithms implemented:

        scaling: Generate a random cell and random distribution of atoms and
                scale the lattice to separate the atoms.

        stretching: Generating a random cell and random distribution of atoms
                    and stretching their bonds until the distance between any
                    two atoms is always greater than the sum of covalent radius.

        :param composition: (pychemia.Composition)
        :param method: (str)
        :param stabilization_number: (int)
        :param nparal: (int)
        :param periodic: (bool)
        :param factor_optimal_volume: (float)
        :return:

        >>> import os
        >>> st = Structure.random_cell('LiAlCl4', stabilization_number=3)
        >>> st.natom
        6
        >>> st.save_json('test.json')
        >>> st2 = Structure.load_json('test.json')
        >>> st == st2
        True
        >>> os.remove('test.json')
        """
        comp = Composition(composition)
        pcm_log.debug('Generating a random structure with composition: ' + str(comp.composition))
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

            ret = pool.map(worker_star, zip(repeat(method), repeat(composition), repeat(periodic), args))

            ngood = 0
            for structure in ret:
                if structure is not None:
                    # print('SH:%d Vol:%10.3f Factor:%10.3f' % (stabilization_history,
                    #                                          structure.volume,
                    #                                          structure.volume / optimal_volume))
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

            if best_volume < factor_optimal_volume * optimal_volume:
                break
            trial += 1

            # log.debug('Trial: %4d  Volume: %7.2f  Optimal Volume: %7.2f  Ratio: %5.2f' %
            #          (trial, best_volume, optimal_volume, best_volume/optimal_volume))

        pool.close()

        if best_structure is not None and periodic:
            # Analysis of the quality for the best structure
            rpos = best_structure.reduced
            for i, j in combinations(range(natom), 2):
                distance = best_structure.lattice.minimal_distance(rpos[i], rpos[j])
                covalent_distance = sum(covalent_radius([symbols[i], symbols[j]]))
                if distance < covalent_distance:
                    pcm_log.debug('Covalent distance: %7.4f  Minimal distance: %7.4f  Difference: %7.3e' %
                                  (covalent_distance, distance, covalent_distance - distance))

        best_structure.canonical_form()
        return best_structure

    @staticmethod
    def random_cluster(composition, method='stretching', stabilization_number=20, nparal=5):

        st = Structure.random_cell(composition=composition, method=method, stabilization_number=stabilization_number,
                                   nparal=nparal, periodic=False)
        return Structure(symbols=st.symbols, positions=st.positions, periodicity=False)

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

        :param cell: A matrix with the 3 unit cell
                     vectors
        :return:
        """
        npcell = np.array(cell)
        if npcell.shape == () or npcell.shape == (1,):
            self.cell = npcell * np.eye(3)
        elif npcell.shape == (3,):
            self.cell = np.diag(npcell)
        else:
            self.cell = np.array(cell).reshape((3, 3))
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
            periodicity: (Boolean) a single value means that the structure has that
                        periodicity all along the 3 directions. Otherwise a list
                        of 3 booleans is required
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

        :param reduced:
        :return:
        """
        self.reduced = np.array(reduced).reshape([-1, 3])

    def sort_sites_using_list(self, sorted_indices):
        sorted_indices = np.array([int(x) for x in sorted_indices])
        self.symbols = list(np.array(self.symbols)[sorted_indices])
        self.positions = self.positions[sorted_indices]
        if self.is_periodic:
            self.reduced = self.reduced[sorted_indices]
        if self.vector_info is not None:
            for vi in self.vector_info:
                if self.vector_info[vi] is not None:
                    self.vector_info[vi] = self.vector_info[vi][sorted_indices]

    def sort_sites(self):

        # First: Sort sites using the distance to the origin
        sorted_indices = np.array([np.linalg.norm(self.positions[i]) for i in range(self.nsites)]).argsort()
        # print sorted_indices
        self.sort_sites_using_list(sorted_indices)

        # Second: Sort again using the atomic number
        if len(self.species) > 1:
            sorted_indices = np.array([atomic_number(x) for x in self.symbols]).argsort()
            self.sort_sites_using_list(sorted_indices)

    def sort_axes(self):
        """
        Sort the lattice vectors in decremental order of their size.
        'a' will be the longest lattice vector
        'c' he shortest
        """

        sorted_indices = self.lattice.lengths.argsort()[::-1]
        self.set_cell(self.cell[sorted_indices])
        self.reduced = self.reduced[:, sorted_indices]

    def align_with_axis(self, axis=0, round_decimals=14):
        lattice = self.lattice
        lattice.align_with_axis(axis=axis, round_decimals=round_decimals)
        self.set_cell(lattice.cell)
        self.reduced2positions()

    def align_with_plane(self, axis=2, round_decimals=14):
        lattice = self.lattice
        lattice.align_with_plane(axis=axis, round_decimals=round_decimals)
        self.set_cell(lattice.cell)
        self.reduced2positions()

    def align_inertia_momenta(self):
        ii = self.inertia_matrix()
        eigval, eigvec = np.linalg.eig(ii)
        eigvec = eigvec.T[eigval.argsort()[::-1]].T
        inveigvec = np.linalg.inv(eigvec)
        self.positions = np.dot(inveigvec, self.positions.T).T

    def canonical_form(self):

        if not self.is_periodic:
            self.relocate_to_cm()
            self.align_inertia_momenta()
        self.sort_sites()
        if self.is_periodic:
            self.sort_axes()
            self.align_with_axis()
            self.align_with_plane()
            self.atoms_in_box()
            self.sort_sites()

    def supercell(self, size):
        """
        Creates a supercell, replicating the positions
        of atoms in the x,y,z directions a number of
        size=(nx,ny,nz) times
        """
        new_natom = np.prod(size) * self.natom
        new_symbols = []
        new_positions = np.zeros((new_natom, 3))
        size = np.array(size).astype(int)

        index = 0
        for i in range(size[0]):
            for j in range(size[1]):
                for k in range(size[2]):
                    for n in range(self.natom):
                        new_symbols.append(self.symbols[n])
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

    @property
    def to_dict(self):

        ret = {'natom': self.natom,
               'symbols': self.symbols,
               'periodicity': self.periodicity,
               'positions': self.positions.tolist(),
               'nspecies': len(self.species),
               'formula': self.formula}
        if self.is_periodic:
            ret['cell'] = self.cell.tolist()
            ret['reduced'] = self.reduced.tolist()
            ret['density'] = self.density
        if self.name is not None:
            ret['name'] = self.name
        if self.comment is not None:
            ret['comment'] = self.comment
        if self.sites != range(self.natom):
            ret['sites'] = list(self.sites)
        if self.occupancies != self.natom * [1.0]:
            ret['occupancies'] = self.occupancies
        # if len(self.vector_info) != 1 or self.vector_info['mag_moments'] is not None:
        #    ret['vector_info'] = self.vector_info
        return ret

    def round(self, decimals=6, pos='reduced'):
        self.set_cell(np.around(self.cell, decimals))
        if pos == 'reduced':
            self.set_reduced(np.around(self.reduced, decimals))
            self.reduced2positions()
        else:
            self.set_positions(np.around(self.positions, decimals))
            self.positions2reduced()

    @staticmethod
    def from_dict(structdict):

        natom = structdict['natom']
        symbols = deep_unicode(structdict['symbols'])
        periodicity = structdict['periodicity']
        positions = np.array(structdict['positions'])

        if 'name' in structdict:
            name = structdict['name']
        else:
            name = None
        if 'comment' in structdict:
            comment = structdict['comment']
        else:
            comment = None
        if 'cell' in structdict:
            cell = np.array(structdict['cell'])
        else:
            cell = None
        if 'reduced' in structdict:
            reduced = np.array(structdict['reduced'])
        else:
            reduced = None
        if 'vector_info' in structdict:
            vector_info = structdict['vector_info']
        else:
            vector_info = None
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
        json.dump(self.to_dict, filep, sort_keys=True, indent=4, separators=(',', ': '))
        filep.close()

    @staticmethod
    def load_json(filename):

        filep = open(filename, 'r')
        structdict = deep_unicode(json.load(filep))
        filep.close()
        return Structure.from_dict(structdict)

    def distance2(self, atom1, atom2):
        assert (isinstance(atom1, int))
        assert (isinstance(atom2, int))
        assert (atom1 < self.natom)
        assert (atom2 < self.natom)

        if self.is_periodic:
            return self.lattice.distance2(self.reduced[atom1], self.reduced[atom2])
        else:
            dm = scipy.spatial.distance_matrix(self.positions, self.positions)
            return dm[atom1, atom2]

    def distance_matrix(self):
        if self.is_periodic:
            dm = np.zeros((self.nsites, self.nsites))
            for i in range(self.nsites - 1):
                for j in range(i + 1, self.nsites):
                    d = self.lattice.distance2(self.reduced[i], self.reduced[j], radius=1E10, limits=[1, 1, 1])
                    # print("%d %d - %d" % (i,j, len(d)))
                    dm[i, j] = min([d[x]['distance'] for x in d])
                    dm[j, i] = dm[i, j]
        else:
            dm = scipy.spatial.distance_matrix(self.positions, self.positions)
        return dm

    def valence_electrons(self):
        ret = 0
        for key, value in self.composition.items():
            ret += value * valence(key)
        return ret

    def __eq__(self, other):
        if self.natom != other.natom:
            ret = False
        elif not np.array_equal(self.positions, other.positions):
            ret = False
        elif not np.array_equal(self.periodicity, other.periodicity):
            ret = False
        elif self.is_periodic and other.is_periodic:
            if not np.array_equal(self.reduced, other.reduced):
                ret = False
            elif not np.array_equal(self.cell, other.cell):
                ret = False
            else:
                ret = True
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
        return self.natom == self.nsites and min(self.occupancies) == 1.0

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
        if self.is_periodic:
            return abs(np.linalg.det(self.cell))
        else:
            volume = (np.max(self.positions[:, 0]) - np.min(self.positions[:, 0])) * \
                     (np.max(self.positions[:, 1]) - np.min(self.positions[:, 1])) * \
                     (np.max(self.positions[:, 2]) - np.min(self.positions[:, 2]))

            if volume > 0.0:
                return volume
            else:
                return 4.0 / 3.0 * np.pi * np.max(self.positions.flatten()) ** 3

    @property
    def species(self):
        return self.get_composition().species

    @property
    def nspecies(self):
        return len(self.get_composition().species)

    @property
    def nsites(self):
        return len(self.positions)

    def scale(self, tolerance=0.7):
        assert self.is_perfect
        assert self.is_crystal
        lattice = self.lattice.scale(self.symbols, self.reduced, tolerance=tolerance)
        return Structure(cell=lattice.cell, reduced=self.reduced, symbols=self.symbols)

    # Not working
    # def cut_void(self, factor=1.5):
    #     ret = self.copy()
    #     ret.canonical_form()
    #     mins = [min(ret.reduced[:, i]) for i in range(3)]
    #     ret.reduced = ret.reduced - mins
    #     ret.reduced2positions()
    #     max_lenght = max(ret.positions[:, 0]) + 2.0*max([covalent_radius(ret.symbols[i]) for i in range(ret.nsites)])
    #     print factor, max_lenght, ret.cell[0,0]
    #     if factor * max_lenght < ret.cell[0, 0]:
    #         cell = ret.cell
    #         cell[0, 0] = factor * max_lenght
    #         ret.set_cell(cell)
    #         return ret
    #     else:
    #         return self

    def atoms_in_box(self):
        while min(self.reduced.flatten()) < 0.0 or max(self.reduced.flatten()) > 1.0:
            self.reduced = (self.reduced + 1.0) % 1.0
            self.reduced2positions()

    def moment_of_inertia(self, axis):

        assert self.is_perfect
        mofi = 0
        for isite in self:
            mofi += mass(isite.symbols[0]) * (sum(np.array(isite.position) ** 2) - isite.position[axis] ** 2)
        return mofi

    def product_of_inertia(self, axis):

        assert self.is_perfect
        pofi = 0
        for isite in self:
            pofi += mass(isite.symbols[0]) * (np.prod(isite.position) / isite.position[axis])
        return pofi

    def inertia_matrix(self):

        im_xx = self.moment_of_inertia(0)
        im_yy = self.moment_of_inertia(1)
        im_zz = self.moment_of_inertia(2)
        im_xy = self.product_of_inertia(2)
        im_xz = self.product_of_inertia(1)
        im_yz = self.product_of_inertia(0)

        im = np.array([[im_xx, -im_xy, -im_xz], [-im_xy, im_yy, -im_yz], [-im_xz, -im_yz, im_zz]])
        return im

    def signature(self):
        comp = self.get_composition()
        gcd = self.get_composition().gcd
        ret = '%02X_%014X_%02X_' % (self.valence_electrons() / gcd, comp.species_hex(), gcd)

        formula = "%s" % comp.sorted_formula(sortby='electroneg')
        formula += (17 - len(formula)) * '_'
        ret += formula

        return ret

    def add_vacuum(self, length, direction=2):

        vacuum = np.zeros(3)
        vacuum[direction] = length
        alpha = self.lattice.alpha
        beta = self.lattice.beta
        gamma = self.lattice.gamma
        newlenghts = self.lattice.lengths + vacuum
        a = newlenghts[0]
        b = newlenghts[1]
        c = newlenghts[2]
        newlattice = self.lattice.from_parameters_to_cell(a, b, c, alpha, beta, gamma)
        return Structure(symbols=self.symbols, cell=newlattice.cell, positions=self.positions)


def load_structure_json(filename):
    ret = Structure()
    ret.load_json(filename)
    return ret


class SiteSet:
    """
    Collection of atomic sites.
    Starting from a Structure the object will create a set of Sites storing the species, occupancies and
    positions of each site.
    """
    def __init__(self, structure):
        """
        SiteSet is a container for a list of Sites

        :param structure: structure from which the Sites will be created.
        """
        self.structure = structure
        self.sitelist = []
        reduced = None

        for isite in range(structure.nsites):
            if structure.sites.count(isite) > 1:
                symbols = []
                occupancies = []
                for jatom in range(structure.natom):
                    if structure.sites[jatom] == isite:
                        symbols.append(structure.symbols[jatom])
                        occupancies.append(structure.occupancies[jatom])
                position = structure.positions[isite]
                if self.structure.is_periodic:
                    reduced = structure.reduced[isite]
            else:
                symbols = [structure.symbols[isite]]
                occupancies = [structure.occupancies[isite]]
                position = structure.positions[isite]
                if self.structure.is_periodic:
                    reduced = structure.reduced[isite]
            self.sitelist.append(Site(symbols=symbols, occupancies=occupancies, position=position))

    def __iter__(self):
        return iter(self.sitelist)


class Site:
    """
    A site is a mapping of one location of space with one or more atoms and the corresponding occupancies.
    The sum of occupancies must be lower or equal to 1.
    """
    def __init__(self, symbols, occupancies=None, position=None):
        """
        Create a atomic site with one or more atoms defined on a unique location on space.

        :param symbols: atomic symbol or list of atomic symbols
        :param occupancies: list of floats with the probability of occupancy for the atom in the site
        :param position: iterable with the location of the atom in space in cartesian coordinates

        .. notes: This class is used for internal operations inside Structure. Not much reason to expose it
        out of this module.

        >>> sit = Site('H')
        >>> sit
        Site(symbols=['H'],occupancies=[1.0],position=[0.0, 0.0, 0.0])
        >>> sit = Site(symbols = ['O', 'N'], occupancies=[0.5, 0.5], position=[1,1,1])
        >>> sit
        Site(symbols=['O', 'N'],occupancies=[0.5, 0.5],position=[1.0, 1.0, 1.0])
        """
        if symbols in atomic_symbols:
            self.symbols = [symbols]
        else:
            try:
                self.symbols = list(symbols)
            except TypeError:
                raise TypeError("%s is not iterable" % symbols)
        for i in self.symbols:
            assert i in atomic_symbols

        if occupancies is None:
            occupancies = 1.0
        if position is None:
            position = 3*[0.0]

        try:
            self.occupancies = [float(x) for x in occupancies]
        except TypeError:
            self.occupancies = [float(occupancies)]

        assert (len(self.occupancies) == len(self.symbols))
        assert sum(self.occupancies) <= 1

        self.position = [float(x) for x in position]

    def __repr__(self):
        ret = self.__class__.__name__+'(symbols=' + repr(self.symbols)
        ret += ',occupancies=' + repr(self.occupancies)
        ret += ',position=' + repr(self.position)
        ret += ')'
        return ret

    def __str__(self):
        return repr(self)


def worker_star(x):
    return random_structure(*x)


def random_structure(method, composition, periodic=True, max_volume=1E10):
    """
    Random Structure created  by random positioning of atoms followed by either scaling of the cell or
    adding a sheer stretching along the smaller distances. The purpose of the lattice change is to avoid any two
    atoms to be closer than the sum of their covalent radius.

    :param method: Can be 'stretching' or 'scaling'.
    :param composition: Can be a Composition object or formula.
    :param periodic: If True, the structure will be periodical in all directions, otherwise a finite system is created.
    :param max_volume: Threshold for creating the Structure, if the volume exceeds the target the method returns None
    :return: Structure if the volume is below than best_volume, None otherwise

    >>> st = random_structure(method='scaling', composition='H2O', periodic=False)
    >>> st.natom
    3
    >>> st = random_structure(method='stretching', composition='NaCl', periodic=True)
    >>> st.natom
    2
    """
    comp = Composition(composition)
    natom = comp.natom
    symbols = comp.symbols
    np.random.seed(struct.unpack("<L", os.urandom(4))[0])

    if periodic:
        new_structure = None

        assert (method in ['scaling', 'stretching'])

	
        while True:
            trial=0
    
            if method == 'scaling':
                lattice = Lattice.random_cell(comp)
                # Random reduced positions
                rpos = np.random.rand(natom, 3)
                mins = [min(rpos[:, i]) for i in range(3)]
                rpos -= mins

                new_lattice = lattice
            else:
                lattice = Lattice.random_cell(comp)
                # Random reduced positions
                rpos = np.random.rand(natom, 3)
                mins = [min(rpos[:, i]) for i in range(3)]
                rpos -= mins

                new_lattice = lattice.stretch(symbols, rpos, tolerance=1.0, extra=0.1)

            if new_lattice.volume < max_volume:
                test = True
                for i in range(natom):
                    for j in range(i + 1, natom):
                        distance = new_lattice.minimal_distance(rpos[i], rpos[j])
                        covalent_dim = sum(covalent_radius([symbols[i], symbols[j]]))
                        if distance < covalent_dim:
                            test = False
                if test:
                    new_structure = Structure(symbols=symbols, reduced=rpos, cell=new_lattice.cell, periodicity=True)
                    minimal_distance = np.min(new_structure.distance_matrix() + 10 * np.eye(new_structure.natom))
                    break
                else:
                    print("Trial failed, distance %f is less than covalent radious %f" % (distance, covalent_dim))
                    trial += 1

            if trial > 100:
                print("Leaving after 100 attemps")
                break

        # else:
        #     print('Volume of Structure %f is larger than max_volume=%f' % (new_lattice.volume, max_volume))
        #     new_structure = None
    else:
        pos = np.random.rand(natom, 3)

        mindis = cluster_minimal_distance(pos)
        if mindis == 0:
            raise ValueError("Distance too small")

        max_cov = np.max(covalent_radius(symbols))
        pos *= max_cov / mindis

        current_volume = (max(pos[:, 0]) - min(pos[:, 0])) * (max(pos[:, 1]) - min(pos[:, 1])) * (
            max(pos[:, 2]) - min(pos[:, 2]))

        if current_volume < max_volume:
            new_structure = Structure(symbols=symbols, positions=pos, periodicity=False)
        else:
            print('Volume of Structure %f is larger than max_volume=%f' % (current_volume, max_volume))
            new_structure = None

    return new_structure
    # End of Worker


def cluster_minimal_distance(pos):
    pos = np.array(pos).reshape((-1, 3))
    dismat = scipy.spatial.distance_matrix(pos, pos)
    tmp = np.max(dismat.flatten())
    return np.min((dismat + tmp * np.eye(len(pos))).flatten())
