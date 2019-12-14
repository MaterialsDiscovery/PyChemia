
import spglib as spg
from pychemia.utils.serializer import generic_serializer
from pychemia.utils.computing import deep_unicode


def spglib_version():
    from pychemia import Structure
    from . import CrystalSymmetry

    # Testing version of spglib
    st = Structure(symbols=['H'])
    symm = CrystalSymmetry(st)
    ret = spg.spglib.spg.dataset(symm.transposed, symm.reduced, symm.numbers, 1e-5, -1.0)
    if type(ret[3]) is list:
        HAS_SPGLIB = False
        version = "%d.%d.%d" % spg.get_version()
        print('SPGLIB current version is %s, please install spglib > 1.9' % version)
    else:
        HAS_SPGLIB = True
    return HAS_SPGLIB


class CrystalSymmetry(object):
    """
    Takes a pychemia.Structure object and creates an object with methods
    to identify symmetry groups and other symmetry related operations.
    Uses spglib to perform various symmetry finding operations.
    """

    def __init__(self, structure):
        """
        Creates a new StructureSymmetry object for a given structure
        This class allows interaction with spglib for several operations related
        to symmetry

        :param structure:

        >>> import pychemia
        >>> a = 4.05
        >>> b = a/2
        >>> fcc = pychemia.Structure(symbols=['Au'], cell=[[0, b, b], [b, 0, b], [b, b, 0]], periodicity=True)
        >>> symm = pychemia.crystal.CrystalSymmetry(fcc)
        >>> symm.number()
        225
        >>> symm.symbol() == u'Fm-3m'
        True

        """
        assert structure.is_crystal
        assert structure.is_perfect
        self.structure = structure
        # Spglib"s convention for the lattice definition is the transpose of cell
        # self._transposed_cell = structure.cell.transpose().copy()
        self.spglib_lattice = generic_serializer(structure.cell.copy())
        # Spglib requires numpy floats.
        # self._transposed_cell = np.array(self._transposed_cell, dtype='double', order='C')
        # self._reduced = np.array(structure.reduced, dtype='double', order='C')
        self.spglib_positions = generic_serializer(structure.reduced)
        # Get a list of indices for each atom in structure
        # Indices starting in 1
        # self._numbers = np.array([structure.species.index(x) + 1 for x in structure.symbols], dtype='intc')
        self.spglib_numbers = [structure.species.index(x) + 1 for x in structure.symbols]

    # @property
    # def transposed(self):
    #     return self._transposed_cell
    #
    # @property
    # def reduced(self):
    #     return self._reduced
    #
    # @property
    # def numbers(self):
    #     return self._numbers

    @property
    def spglib_cell(self):
        return self.spglib_lattice, self.spglib_positions, self.spglib_numbers

    def get_spacegroup(self, symprec=1e-5):
        return spg.get_spacegroup(cell=self.spglib_cell, symprec=symprec)

    def get_symmetry(self, symprec=1e-5):
        return spg.get_symmetry(cell=self.spglib_cell, symprec=symprec)

    def get_symmetry_dataset(self, symprec=1e-5):
        ret = spg.get_symmetry_dataset(cell=self.spglib_cell, symprec=symprec)
        if ret is None:
            raise ValueError(self.spglib_cell)
        return ret

    def number(self, symprec=1e-5):
        return int(self.get_symmetry_dataset(symprec)['number'])

    def symbol(self, symprec=1e-5):
        return deep_unicode(self.get_symmetry_dataset(symprec)['international'])

    def hall_number(self, symprec=1e-5):
        return self.get_symmetry_dataset(symprec)['hall_number']

    def get_new_structure(self, spglib_cell):
        from pychemia import Structure
        cell = spglib_cell[0]
        reduced = spglib_cell[1]
        symbols = [self.structure.species[x - 1] for x in spglib_cell[2]]
        return Structure(cell=cell, reduced=reduced, symbols=symbols)

    def refine_cell(self, symprec=1e-5):
        """
        Refine a pychemia Structure using the tolerances and return a new structure in a Bravais lattice

        :param symprec: (float) Tolerance of distance between atomic positions and between lengths of lattice vectors
                        to be tolerated in the symmetry finding.
        :return: A new pychemia Structure in a Bravais lattice

        :rtype : (pychemia.Structure)

        Example:
        >>> import pychemia
        >>> a = 4.05
        >>> b = a/2
        >>> fcc = pychemia.Structure(symbols=['Au'],
        ...       cell=[[0, b+1E-5, b-1E-5], [b+1E-5, 0, b-1E-5], [b+1E-5, b-1E-5, 0]], periodicity=True)
        >>> symm = pychemia.crystal.CrystalSymmetry(fcc)
        >>> symm.number()
        2
        >>> symm.symbol() == u'P-1'
        True
        >>> fcc2 = symm.refine_cell(symprec=1E-3)
        >>> symm2 = pychemia.crystal.CrystalSymmetry(fcc2)
        >>> symm2.number()
        225
        >>> symm2.symbol() == u'Fm-3m'
        True

        """
        new_spglib_cell = spg.refine_cell(self.spglib_cell, symprec=symprec)
        return self.get_new_structure(new_spglib_cell)

    def find_primitive(self, symprec=1e-5):
        new_spglib_cell = spg.find_primitive(self.spglib_cell, symprec=symprec)
        return self.get_new_structure(new_spglib_cell)

    def standarize_cell(self, symprec=1e-5):
        new_spglib_cell = spg.standarize_cell(self.spglib_cell, symprec=symprec)
        return self.get_new_structure(new_spglib_cell)

    def get_space_group_type(self, symprec=1e-5):
        return spg.get_spacegroup_type(self.hall_number(symprec=symprec))

    def crystal_system(self, symprec=1e-5):

        num = self.number(symprec)

        if num < 3:
            return u'Triclinic'
        elif num < 16:
            return u'Monoclinic'
        elif num < 75:
            return u'Orthorhombic'
        elif num < 143:
            return u'Tetragonal'
        elif num < 168:
            return u'Trigonal'
        elif num < 195:
            return u'Hexagonal'
        else:
            return u'Cubic'

    def symmetrize(self, initial_symprec=0.01, final_symprec=0.1, delta_symprec=0.01):
        if self.structure.natom == 1:
            return self.structure.copy()
        prec = initial_symprec
        while prec < final_symprec:
            if self.number(symprec=prec) > self.number(symprec=initial_symprec):
                break
            else:
                prec += delta_symprec
        if prec > final_symprec:
            prec = final_symprec
        new_bravais = self.refine_cell(symprec=prec)
        sym2 = CrystalSymmetry(new_bravais)
        return sym2.find_primitive()

        # def get_symmetry_dataset(self, symprec=1e-5, angle_tolerance=-1.0):
        #
        #     keys = ('number',
        #             'hall_number',
        #             'international',
        #             'hall',
        #             'transformation_matrix',
        #             'origin_shift',
        #             'rotations',
        #             'translations',
        #             'wyckoffs',
        #             'equivalent_atoms',
        #             'brv_lattice',
        #             'brv_types',
        #             'brv_positions')
        #     dataset = {}
        #     for key, data in zip(keys, spg.spglib.spg.dataset(self._transposed_cell,
        #                                                       self._reduced,
        #                                                       self._numbers,
        #                                                       symprec,
        #                                                       angle_tolerance)):
        #         dataset[key] = data
        #
        #     dataset['international'] = dataset['international'].strip()
        #     dataset['hall'] = dataset['hall'].strip()
        #     dataset['transformation_matrix'] = np.array(
        #         dataset['transformation_matrix'], dtype='double', order='C')
        #     dataset['origin_shift'] = np.array(dataset['origin_shift'], dtype='double')
        #     dataset['rotations'] = np.array(dataset['rotations'],
        #                                     dtype='intc', order='C')
        #     dataset['translations'] = np.array(dataset['translations'],
        #                                        dtype='double', order='C')
        #     letters = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ"
        #     dataset['wyckoffs'] = [letters[x] for x in dataset['wyckoffs']]
        #     dataset['equivalent_atoms'] = np.array(dataset['equivalent_atoms'],
        #                                            dtype='intc')
        #     dataset['brv_lattice'] = np.array(np.transpose(dataset['brv_lattice']),
        #                                       dtype='double', order='C')
        #     dataset['brv_types'] = np.array(dataset['brv_types'], dtype='intc')
        #     dataset['brv_positions'] = np.array(dataset['brv_positions'],
        #                                         dtype='double', order='C')
        #
        #     return dataset
        #
        # def get_spacegroup(self, symprec=1e-5, angle_tolerance=-1.0):
        #     """
        #     Return space group in international table symbol and number
        #     as a string.
        #
        #     :param symprec: (float) Symmetry precision
        #     :param angle_tolerance: (float) angle tolerance for spglib internal routine
        #     :return:
        #     """
        #
        #     dataset = self.get_symmetry_dataset(symprec=symprec, angle_tolerance=angle_tolerance)
        #     return "%s (%d)" % (dataset['international'], dataset['number'])
        #
        # def spacegroup(self, symprec=1e-5, angle_tolerance=-1.0):
        #     """
        #     Computes the space group for the structure with a given
        #     precision in distances (symprec) and angle tolerance in degrees (angle_tolerance)
        #
        #     :param symprec: (float) Tolerance of distance between atomic positions and between lengths of
        #                             lattice vectors
        #     to be tolerated in the symmetry finding.
        #     :param angle_tolerance: (float) Tolerance of angle between lattice vectors in degrees to be
        #                                     tolerated in the symmetry finding.
        #
        #     :return: The space group with the symbol and number as a string
        #     :rtype : (str)
        #     """
        #
        #     return self.get_spacegroup(symprec, angle_tolerance)
        #
        # def symbol(self, symprec=1e-5, angle_tolerance=-1.0):
        #     """
        #     Computes the space group symbol for the structure with a given
        #     precision in distances (symprec) and angle tolerance in degrees (angle_tolerance)
        #
        #     :param symprec: (float) Tolerance of distance between atomic positions and between lengths of
        #                             lattice vectors to be tolerated in the symmetry finding.
        #     :param angle_tolerance: (float) Tolerance of angle between lattice vectors in degrees to be
        #                                     tolerated in the symmetry finding.
        #
        #     :return: The space group symbol as a string
        #     :rtype : (str)
        #     """
        #     return self.spacegroup(symprec, angle_tolerance).split()[0]
        #
        # def number(self, symprec=1e-5, angle_tolerance=-1.0):
        #     """
        #     Computes the space group number for the structure with a given
        #     precision in distances (symprec) and angle tolerance in degrees (angle_tolerance)
        #
        #     :param symprec: (float) Tolerance of distance between atomic positions and between lengths
        #                             of lattice vectors to be tolerated in the symmetry finding.
        #     :param angle_tolerance: (float) Tolerance of angle between lattice vectors in degrees to be
        #                             tolerated in the symmetry finding.
        #
        #     :return: The space group number
        #     :rtype : (int)
        #     """
        #     return int(self.spacegroup(symprec, angle_tolerance).split()[-1].strip()[1:-1])

        # def refine_cell(self, symprec=1e-5, angle_tolerance=-1.0):
        #     from pychemia import Structure
        #     natom = self.structure.natom
        #     cell = np.array(self.structure.cell.transpose(), dtype='double', order='C')
        #
        #     pos = np.zeros((natom * 4, 3), dtype='double')
        #     pos[:natom] = self._reduced
        #
        #     numbers = np.zeros(natom * 4, dtype='intc')
        #     numbers[:natom] = np.array(self._numbers, dtype='intc')
        #
        #     natom_bravais = spg.spglib.spg.refine_cell(cell, pos, numbers, natom, symprec, angle_tolerance)
        #     if natom_bravais == 0:
        #         return self.structure.copy()
        #
        #     cell = np.array(cell.T, dtype='double', order='C')
        #     reduced = np.array(pos[:natom_bravais], dtype='double', order='C')
        #     symbols = [self.structure.species[x] for x in (numbers[:natom_bravais] - 1)]
        #
        #     return Structure(cell=cell, symbols=symbols, reduced=reduced)
        #
        # def find_primitive(self, symprec=1e-5, angle_tolerance=1.0):
        #     """
        #     Search the primitive pychemia Structure using the tolerances.
        #     If no primitive cell is found a copy of the original structure is returned
        #
        #     :param symprec: (float) Tolerance of distance between atomic positions and between lengths of
        #                             lattice vectors to be tolerated in the symmetry finding.
        #     :param angle_tolerance: (float) Tolerance of angle between lattice vectors in degrees to
        #                                     be tolerated in the symmetry finding.
        #
        #     :return: A new pychemia Structure in a Bravais lattice
        #     :rtype : (pychemia.Structure)
        #     """
        #     from pychemia import Structure
        #     # Create copies of the arguments
        #     cell = np.array(self._transposed_cell, dtype='double', order='C')
        #     reduced = np.array(self._reduced, dtype='double', order='C')
        #     numbers = np.array(self._numbers, dtype='intc')
        #
        #     natom_prim = spg.spglib.spg.primitive(cell, reduced, numbers, symprec, angle_tolerance)
        #     symbols = [self.structure.species[x] for x in (numbers[:natom_prim] - 1)]
        #     reduced = reduced[:natom_prim]
        #
        #     if natom_prim > 0:
        #         return Structure(cell=cell.T, reduced=reduced, symbols=symbols)
        #     else:
        #         return self.structure.copy()
        #
