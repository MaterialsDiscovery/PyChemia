"""
Definition of the class Structure
This class defines methods to create and manipulate
atomic structures such as molecules, clusters and crystals
"""
import itertools

import numpy as _np
import json as _json
from math import sin, cos, exp, sqrt
import math
import sys
from pychemia.geometry.lattice import Lattice
from pychemia.geometry.delaunay import get_reduced_bases
from pychemia.geometry.composition import Composition
from pychemia.utils.computing import unicode2string
from pychemia.utils.mathematics import distances
from pychemia.utils.periodic import mass, atomic_number, covalent_radius, valence, atomic_symbol, atomic_symbols

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
    """

    def __init__(self, **kwargs):
        """
        Structure

        Represents a molecule, cluster, wire, slab or crystal structure
        The positions of the atoms and their atomic symbols are declared
        in 'positions' and 'symbols' respectively.

        For periodic structures, the 'periodicity' can be declared.
        and cell parameters in 'cell'

        Magnetic moments can be associated in the array vector_info['magnetic_moments'].

        This object contains no dynamical information. That information
        is supported by the child class DynamicStructure

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
        symmetry   : Symmetry information like point symmetry operations
                     and space group
        name       : Free text to identify structure (Only one line, max 50 chars)
        comment    : Free text to identify structure

        Examples:

        >>> import pychemia
        >>> a = pychemia.geometry.Structure()
        >>> print a
        Empty structure
        >>> a = pychemia.geometry.Structure(symbols=['Xe'])
        >>> print a.natom
        1
        >>> d = 1.104
        >>> a = pychemia.geometry.Structure(symbols=['N', 'N'], positions=[[0, 0, 0], [0, 0, d]], periodicity=False)
        >>> print a.natom
        2
        >>> a = 4.05
        >>> b = a/2
        >>> fcc = pychemia.geometry.Structure(symbols=['Au'], cell=[[0, b, b], [b, 0, b], [b, b, 0]], periodicity=True)
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

        self._lattice = None
        self._composition = None

        # Fill the values from args
        if 'name' in kwargs:
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
            cell = _np.array(kwargs['cell'])
            self.set_cell(cell)
        if 'positions' in kwargs:
            positions = _np.array(kwargs['positions'])
            self.set_positions(positions)
        if 'reduced' in kwargs:
            reduced = _np.array(kwargs['reduced'])
            self.set_reduced(reduced)
        if 'mag_moments' in kwargs:
            self.set_mag_moments(_np.array(kwargs['mag_moments']))

        # Lets autocomplete the missing information
        self._autocomplete()

        if not self._check():
            print('Arguments non consistent')

    def __len__(self):
        """
        Return the number of atoms in the cell
        """
        return self.natom

    def __str__(self):
        """
        String representation of the object
        """
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
        ret = 'Structure(symbols='+str(self.symbols)
        if self.is_periodic:
            if _np.all(_np.diag(self.cell.diagonal()) == self.cell):
                if _np.max(self.cell.diagonal()) == _np.min(self.cell.diagonal()):
                    ret += ', cell='+str(self.cell[0, 0])
                else:
                    ret += ', cell='+str(self.cell.diagonal().tolist())
            else:
                ret += ', cell='+str(self.cell.tolist())
        ret += ', positions='+str(self.positions.tolist())
        if all([self.periodicity[0] == item for item in self.periodicity]):
            ret += ', periodicity='+str(self.periodicity[0])
        else:
            ret += ', periodicity='+str(self.periodicity)
        ret += ')'
        return ret

    @property
    def is_periodic(self):
        return any(self.periodicity)

    def _autocomplete(self):
        """
        Autocomplete items in the structure
        that could be obtain from other information
        present in the object
        """
        if self.natom is None:
            if not self.positions is None:
                self.natom = len(self.positions)
            elif not self.reduced is None:
                self.natom = len(self.reduced)
            elif not self.symbols is None:
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
            if not self.reduced is None:
                self.reduced2positions()
            else:
                if self.natom == 0:
                    self.positions = _np.array([])
                elif self.natom == 1:
                    self.positions = _np.array([[0.0, 0.0, 0.0]])
                else:
                    raise ValueError('Positions must be present for more than 1 atom')

        if self.reduced is None and self.is_crystal:
            if self.positions is not None and self.natom > 0:
                self.positions2reduced()
            else:
                self.reduced = _np.array([])

    @property
    def is_crystal(self):
        if not self.is_periodic:
            return False
        else:
            return self.get_cell().periodic_dimensions == 3

    def _check(self):
        """
        Internal check of the selfconsistency of the
        atomic structure
        """
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
        assert(name in atomic_symbols)
        assert(option in ['cartesian', 'reduced'])
        self.symbols = self.symbols.append(name)
        self.natom += 1

        if option == 'cartesian':
            if self.natom == 0:
                self.positions = _np.array(coordinates).reshape([-1, 3])
            else:
                self.positions = _np.append(self.positions, coordinates).reshape([-1, 3])
            self.positions2reduced()
        elif option == 'reduced':
            if self.natom == 0:
                self.reduced = _np.array(coordinates).reshape([-1, 3])
            else:
                self.reduced = _np.append(self.reduced, coordinates).reshape([-1, 3])
            self.reduced2positions()

    def del_atom(self, index):
        """
        Removes the atom with the given index

        :param index:
        :return:
        """
        assert(abs(index) < self.natom)
        self.symbols.pop[index]
        _np.delete(self.positions, index, 0)
        _np.delete(self.reduced, index, 0)
        self.natom -= 1

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
        center_of_mass = _np.zeros(3)
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

        rotationx = _np.array([[1, 0, 0], [0, cos(tx), -sin(tx)], [0, sin(tx), cos(tx)]])
        rotationy = _np.array([[cos(ty), 0, sin(ty)], [0, 1, 0], [-sin(ty), 0, cos(ty)]])
        rotationz = _np.array([[cos(tz), -sin(tz), 0], [sin(tz), cos(tz), 0], [0, 0, 1]])

        rotation = _np.dot(_np.dot(rotationx, rotationy), rotationz)

        for i in range(self.natom):
            self.positions[i] = _np.dot(rotation, self.positions[i])

    @property
    def composition(self):
        return self.get_composition().composition

    @property
    def formula(self):
        return self.get_composition().formula

    def get_cell(self):
        if self._lattice is None:
            self._lattice = Lattice(self.cell)
        return self._lattice

    def get_composition(self, gcd=True):
        """
        Computes the composition of the Structure
        as the count of each species in the cell
        If gcd is True the values are divided by the 
        greatest common divisor 

        :param gcd: bool

        :rtype : dict
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

    @property
    def density(self):
        """
        Computes the density of the cell

        :rtype : float
        """
        return sum(_np.array(mass(self.symbols))) / self.volume

    @property
    def volume(self):
        """
        Computes the volume of the cell

        :rtype : float
        """
        return abs(_np.linalg.det(self.cell))

    def positions2reduced(self):
        """
        Computes the cell-reduced coordinates from the
        cartesian dimensional coordinates
        """
        self.reduced = _np.linalg.solve(self.cell.T, self.positions.T).T
        for i in range(3):
            if self.periodicity[i]:
                self.reduced[:, i] %= 1.0

    def reduced2positions(self):
        """
        Computes the dimensional cartesian coordinates
        from the adimensional cell-reduced coordinates
        """
        self.positions = _np.dot(self.reduced, self.cell)

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

        Args:
            iatom: (int) index of first atom
            jatom: (int) index of second atom
            with_periodicity: (bool) if the periodic
                              images should be considered to compute
                              the shortest distance
            tolerance: (float) Tolerance for the bases reduction
        """

        if with_periodicity:
            reduced_bases = get_reduced_bases(self.cell, tolerance)
            scaled_pos = _np.dot(self.positions, _np.linalg.inv(reduced_bases))
            # move scaled atomic positions into -0.5 < r <= 0.5
            for pos in scaled_pos:
                pos -= pos.round()

            # Look for the shortest one in surrounded 3x3x3 cells
            distances_list = []
            for i in (-1, 0, 1):
                for j in (-1, 0, 1):
                    for k in (-1, 0, 1):
                        distances_list.append(_np.linalg.norm(
                            _np.dot(scaled_pos[iatom] - scaled_pos[jatom] +
                                    _np.array([i, j, k]), reduced_bases)))
            ret = min(distances_list)

        else:
            posi = self.positions[iatom]
            posj = self.positions[jatom]
            ret = _np.linalg.norm(posi - posj)

        return ret

    def set_cell(self, cell):
        """
        Set the vectors defining the cell

        Args:
            cell: A matrix with the 3 unit cell
            vectors
        """
        npcell = _np.array(cell)
        if npcell.shape == () or npcell.shape == (1,):
            self.cell = npcell*_np.eye(3)
        elif npcell.shape == (3,):
            self.cell = _np.diag(npcell)
        else:
            self.cell = _np.array(cell).reshape([3, 3])

    def set_mag_moments(self, mag_moments):
        """
        Set the magnetic moments with one vector on each
        atom

        Args:
            mag_moments: List or numpy array with one
            vector for each atom.
            The values will be converted into a numpy array
        """
        self.vector_info['mag_moments'] = _np.array(mag_moments).reshape([-1, 3])

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
        self.positions = _np.array(positions).reshape([-1, 3])

    def set_reduced(self, reduced):
        """
        Set the reduced positions of the atoms
        This contains adimensional values
        relative to cell vectors

        Args:
            positions: A array of 3 vectors
            with adimensional coordinates
        """
        self.reduced = _np.array(reduced).reshape([-1, 3])

    def sort_byaxis(self, axis):
        """
        Sort the atoms by the given axis
        """
        if axis == 'x':
            index = 0
        elif axis == 'y':
            index = 1
        elif axis == 'z':
            index = 2
        else:
            return
        order = _np.argsort(self.positions[:, index])
        self.positions = self.positions[order]
        self.symbols = self.symbols[order]

    def supercell(self, nx, ny, nz):
        """
        Creates a supercell, replicating the positions
        of atoms in the x,y,z directions a number of
        nx,ny,nz times
        """
        old_natom = self.natom
        for n in range(old_natom):
            for i in range(nx):
                for j in range(ny):
                    for k in range(nz):
                        if (i + j + k) > 0:
                            name = self.symbols[n]
                            position = self.positions[n] + (i * self.cell[0] + j * self.cell[1] + k * self.cell[2])
                            self.add_atom(name, position)
        self.cell[0] = nx * self.cell[0]
        self.cell[1] = ny * self.cell[1]
        self.cell[2] = nz * self.cell[2]

    def copy(self, deep=False):
        """
        Get a copy of the object
        """
        copy_struct = Structure(cell=self.cell,
                                positions=self.positions,
                                periodicity=self.periodicity,
                                symbols=self.symbols)
        return copy_struct

    def get_all_distances(self):
        bonds_dict = {}
        index = 0
        all_distances = []
        for i, j in itertools.combinations(range(self.natom), 2):
            ret = self.get_cell().distance2(self.reduced[i], self.reduced[j])
            for k in ret:
                if str(i) not in bonds_dict:
                    bonds_dict[str(i)] = [index]
                else:
                    bonds_dict[str(i)].append(index)
                if str(j) not in bonds_dict:
                    bonds_dict[str(j)] = [index]
                else:
                    bonds_dict[str(j)].append(index)
                ret[k]['pair'] = (i, j)
                all_distances.append(ret[k])
                index += 1
        for i in range(self.natom):
            ret = self.get_cell().distance2(self.reduced[i], self.reduced[i])
            for k in ret:
                if str(i) not in bonds_dict:
                    bonds_dict[str(i)] = [index]
                else:
                    bonds_dict[str(i)].append(index)
                ret[k]['pair'] = (i, i)
                all_distances.append(ret[k])
                index += 1
        return bonds_dict, all_distances

    def get_bonds_coordination(self, tolerance=1, ensure_conectivity=False):

        bonds_dict, all_distances = self.get_all_distances()
        bonds = []
        tolerances=[]
        for i in range(self.natom):
            tole = tolerance
            while True:
                tmp_bonds = []
                min_proportion = sys.float_info.max
                for j in bonds_dict[str(i)]:
                    atom1 = self.symbols[all_distances[j]['pair'][0]]
                    atom2 = self.symbols[all_distances[j]['pair'][1]]
                    sum_covalent_radius = sum(covalent_radius([atom1, atom2]))
                    distance = all_distances[j]['distance']
                    if distance == 0.0:
                        continue
                    proportion = distance/sum_covalent_radius
                    min_proportion = min(min_proportion, proportion)
                    if proportion <= tole:
                        #print all_distances[j]
                        tmp_bonds.append(j)
                if len(tmp_bonds) == 0 and ensure_conectivity:
                    #print 'Changing tolerance'
                    tole = min_proportion
                else:
                    bonds.append(tmp_bonds)
                    tolerances.append(min_proportion)
                    break

        coordination = [len(x) for x in bonds]
        return bonds, coordination, all_distances, tolerances

    def hardness(self, noupdate=False, verbose=False, tolerance=1):
        """
        Calculates the hardness of a structure based in the model of XX
        We use the covalent radii from pychemia.utils.periodic.
        If noupdate=False
        the Laplacian matrix method is not used and rcut is 2*max(cov_radii)

        :param noupdate: (bool) If True, the Laplacian method is used
        :param verbose: (bool) To print some debug info
        :param tolerance: (float)

        :rtype : (float)
        """

        bonds, coordination, all_distances, tolerances = self.get_bonds_coordination(tolerance=tolerance, ensure_conectivity=True)

        if verbose:
            print 'BONDS'
            print bonds
            print 'COORDINATION'
            print coordination

        sigma = 3.0
        c_hard = 1300.0
        x = 1.
        tot = 0.0
        f_d = 0.0
        f_n = 1.0
        atomicnumbers = atomic_number(self.get_composition().species)

        if verbose:
            print atomicnumbers

        for i in atomicnumbers:
            f_d += valence(i) / covalent_radius(i)
            f_n *= valence(i) / covalent_radius(i)
        f = 1.0 - (len(atomicnumbers) * f_n ** (1.0 / len(atomicnumbers)) / f_d) ** 2

        # Selection of different bonds
        diff_bonds = _np.unique(_np.array(reduce(lambda x, y: x+y, bonds)))
        for i in diff_bonds:
            i1 = all_distances[i]['pair'][0]
            i2 = all_distances[i]['pair'][1]

            ei = valence(self.symbols[i1]) / covalent_radius(self.symbols[i1])
            ej = valence(self.symbols[i2]) / covalent_radius(self.symbols[i2])
            #print 'bond ->', sqrt(ei * ej), (coordination[i1] * coordination[i2]), all_distances[i]['distance']

            sij = sqrt(ei * ej) / (coordination[i1] * coordination[i2]) / all_distances[i]['distance']
            num_i_j_bonds = len([j for j in diff_bonds if i1 in all_distances[j]['pair'] and i2 in all_distances[j]['pair']])
            #print 'sij', sij
            #print 'num_i_j_bonds', num_i_j_bonds
            tot += num_i_j_bonds
            x *= sij
            #print 'x', x

        if verbose:
            print("V:", self.volume)
            print("f:", f)
            print("x:", x)

        #print 'len_bonds', len(diff_bonds)
        #print 'hardness_value =', c_hard, self.volume, (len(diff_bonds)), ( x ** (1. / (len(diff_bonds)))), exp(-sigma * f)
        hardness_value = c_hard / self.volume * (len(diff_bonds) * x ** (1. / (len(diff_bonds)))) * exp(-sigma * f)

        if verbose:
            print hardness_value

        return round(hardness_value, 3)

    def get_bonds(self, radius, noupdate=False, verbose=False, tolerance=0.05):
        """
        Calculates bond lengths between all atoms within "radius".

        Args:
            radius: (float) The radius of the sphere were the bonds
                    will be computed
            noupdate: (bool) If the original radius should be increased
                      until a connected structure is obtained
            verbose: (bool) Print some info about the number of bonds
                     computed
            tolerance: (float) Tolerance for creation of a bond

        Return:
            The values of the bonds are stored in self.info['bonds']

        :param radius: (float)
        :param noupdate: (bool)
        :param verbose: (bool)
        :param tolerance: (float)

        :rtype : dict
        """
        # The number of atoms in the cell
        n = self.natom

        if verbose:
            print('Testing with %s of atoms in cell' % n)
            print('Starting with R_cut = %s' % radius)

        while True:
            lap_m = _np.zeros((n, n))
            dis_dic = {}
            ndifbonds = 0
            found = False
            # lap_mat[i,j] = lap_mat[j,i]
            for i in range(n - 1):
                for j in range(i + 1, n):
                    dis = self.get_distance(i, j, with_periodicity=True)
                    if dis < radius:
                        if len(dis_dic) != 0:
                            for kstr, kj in dis_dic.items():
                                tstr = "%s%s" % (self.symbols[i],
                                                 self.symbols[j])
                                if abs(kj[0] - dis) < tolerance and tstr in kstr:
                                    dis_dic[kstr][1] += 1
                                    found = True
                                    break
                        if not found:
                            ndifbonds += 1
                            kstr = "%s%s%s%s" % (self.symbols[i],
                                                 self.symbols[j],
                                                 self.symbols[i],
                                                 ndifbonds)
                            dis_dic[kstr] = [dis, 1, [i, j]]

                        lap_m[i][j] = -1
                        lap_m[j][i] = -1

            # lap_map[i,i]
            coordination = []
            for i in range(n):
                lap_m[i, i] = abs(sum(lap_m[i]))
                coordination.append(lap_m[i, i])

            # Get eigenvalues and check how many zeros;
            eigen = _np.linalg.eig(lap_m)[0]

            nzeros = 0
            for i in eigen:
                if round(i.real, 5) == 0.0:
                    nzeros += 1
            if nzeros == 1 or (noupdate and len(dis_dic) > 0):
                break
            else:
                radius += 0.02

            if radius > 10.0:
                print("Cut off radius > 10.0")
                break

        if verbose:
            print('New cut off = %s' % radius)
            print('Bonds:')
            for i, j in dis_dic.items():
                print("  %s %s" % (i, j))

        return radius, coordination, dis_dic

    def hardness_OLD(self, noupdate=False, verbose=False, tolerance=0.05):
        """
        Calculates the hardness of a structure based in the model of XX
        We use the covalent radii from pychemia.utils.periodic.
        If noupdate=False
        the Laplacian matrix method is not used and rcut is 2*max(cov_radii)

        :param noupdate: (bool) If True, the Laplacian method is used
        :param verbose: (bool) To print some debug info
        :param tolerance: (float)

        :rtype : (float)
        """

        #from ase.data import covalent_radii

        #atms=atms.repeat([2,2,2])
        spc = self.copy()
        spc.supercell(2, 2, 2)
        natom = spc.natom
        volume = spc.volume

        #max_covalent_radius = max([covalent_radii[i.number] for i in atms])
        max_covalent_radius = max(covalent_radius(spc.symbols))
        if verbose:
            print('Number of atoms', natom)
            print('Volume         ', volume)
            print('Covalent rad max', max_covalent_radius)
        #rcut, coord, dis_dic  = get_bonds(atms,2.0*max_covalent_radius, noupdate,verbose,tolerance)
        rcut, coord, dis_dic = spc.get_bonds(2.0 * max_covalent_radius, noupdate, verbose, tolerance)

        sigma = 3.0
        c_hard = 1300.0
        x = 1.
        tot = 0.0
        f_d = 0.0
        f_n = 1.0
        dic_atms = {}
        #for i in atms:
        #    dic_atms[i.symbol] = i.number
        for i in spc.symbols:
            dic_atms[i] = atomic_number(i)

        #for i in dic_atms.keys():
        #    f_d += Z[i] / covalent_radii[dic_atms[i]]
        #    f_n *= Z[i] / covalent_radii[dic_atms[i]]
        #f = 1.0 - (len(dic_atms)*f_n**(1.0/len(dic_atms)) / f_d)**2

        for i in dic_atms.keys():
            f_d += valence(i) / covalent_radius(i)
            f_n *= valence(i) / covalent_radius(i)
        f = 1.0 - (len(dic_atms) * f_n ** (1.0 / len(dic_atms)) / f_d) ** 2

        if verbose:
            print 'BONDS'
            print dis_dic
            print 'COORDINATION'
            print coord

        for i in dis_dic.keys():
            i1 = dis_dic[i][2][0]
            i2 = dis_dic[i][2][1]

            #ei = Z[atms[i1].symbol] / covalent_radii[atms[i1].number]
            #ej = Z[atms[i2].symbol] / covalent_radii[atms[i2].number]
            ei = valence(spc.symbols[i1]) / covalent_radius(spc.symbols[i1])
            ej = valence(spc.symbols[i2]) / covalent_radius(spc.symbols[i2])


            #print 'bond ->', sqrt(ei * ej), (coord[i1] * coord[i2]), dis_dic[i][0]

            #        print atms[i1].symbol, ei, atms[i2].symbol, ej
            sij = sqrt(ei * ej) / (coord[i1] * coord[i2]) / dis_dic[i][0]
            #print 'sij', sij
            #print 'num_i_j_bonds', dis_dic[i][1]

            tot += dis_dic[i][1]
            x *= sij * dis_dic[i][1]
            #print 'x', x

        if verbose:
            print("V:", volume)
            print("f:", f)
            print("x:", x)

        #print 'len_bonds', len(dis_dic)
        #print 'hardness_value =', c_hard, volume, (len(dis_dic)), (  x ** (1. / (len(dis_dic)))), exp(-sigma * f)
        hardness_value = c_hard / volume * (len(dis_dic) * x ** (1. / (len(dis_dic)))) * exp(-sigma * f)

        if verbose:
            print hardness_value

        return round(hardness_value, 3)

    def plot(self):
        from mayavi import mlab
        assert(self.natom > 0)

        x = self.positions[:, 0]
        y = self.positions[:, 1]
        z = self.positions[:, 2]
        cr = covalent_radius(self.symbols)

        mlab.points3d(x, y, z, cr, scale_factor=1)

        if self.is_crystal:
            frame, line1, line2, line3 = self.get_cell().get_path()

            mlab.plot3d(frame[:, 0], frame[:, 1], frame[:, 2], tube_radius=.05, color=(1, 1, 1))
            mlab.plot3d(line1[:, 0], line1[:, 1], line1[:, 2], tube_radius=.05, color=(1, 1, 1))
            mlab.plot3d(line2[:, 0], line2[:, 1], line2[:, 2], tube_radius=.05, color=(1, 1, 1))
            mlab.plot3d(line3[:, 0], line3[:, 1], line3[:, 2], tube_radius=.05, color=(1, 1, 1))

        mlab.view()
        return mlab.gcf()

    def todict(self):

        ret = {'name': self.name,
               'comment': self.comment,
               'natom': self.natom,
               'symbols': self.symbols,
               'periodicity': self.periodicity,
               'cell': self.cell.tolist(),
               'positions': self.positions.tolist(),
               'reduced': self.reduced.tolist(),
               'vector_info': self.vector_info
               }
        return ret

    def fromdict(self, structdict):

        self.name = structdict['name']
        self.comment = structdict['comment']
        self.natom = structdict['natom']
        self.symbols = structdict['symbols']
        self.periodicity = structdict['periodicity']
        self.cell = _np.array(structdict['cell'])
        self.positions = _np.array(structdict['positions'])
        self.reduced = _np.array(structdict['reduced'])
        self.vector_info = structdict['vector_info']

    def save_json(self, filename):

        filep = open(filename, 'w')
        _json.dump(self.todict(), filep, sort_keys=True, indent=4, separators=(',', ': '))
        filep.close()

    def load_json(self, filename):

        filep = open(filename, 'r')
        structdict = unicode2string(_json.load(filep))
        self.fromdict(structdict)

    def distance2(self, atom1, atom2):
        assert (isinstance(atom1, int))
        assert (isinstance(atom2, int))
        assert (atom1 < self.natom)
        assert (atom2 < self.natom)
        return self.get_cell().distance2(self.reduced[atom1], self.reduced[atom2])

    def __eq__(self, other):
        if self.natom != other.natom:
            ret = False
        elif not _np.array_equal(self.cell, other.cell):
            ret = False
        elif not _np.array_equal(self.positions, other.positions):
            ret = False
        elif not _np.array_equal(self.reduced, other.reduced):
            ret = False
        elif not _np.array_equal(self.periodicity, other.periodicity):
            ret = False
        else:
            ret = True
        return ret

    def __ne__(self, other):
        return not self.__eq__(other)


def load_structure_json(filename):
    ret = Structure()
    ret.load_json(filename)
    return ret


class DynamicStructure(Structure):
    """
    A DynamicStructure contains extra information such
    as velocities, Constrains in the movement of atoms, etc
    """

    def __init__(self, **kwargs):
        Structure.__init__(self, **kwargs)

class MetaStructure():
    """
    This class offer the possibility for atoms of being in a certain
    position with a probability lower than 1
    For example alloys and structural vacancies
    """
    def __init__(self):
        pass
