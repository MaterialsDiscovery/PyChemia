import numpy as np

from pychemia import Structure

try:
    import spglib as spg

    USE_SPGLIB = True
except ImportError:
    USE_SPGLIB = False
    # from pyspglib import spglib as spg


class StructureSymmetry(object):
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

        Example:
>>> import pychemia
>>> a = 4.05
>>> b = a/2
>>> fcc = pychemia.Structure(symbols=['Au'], cell=[[0, b, b], [b, 0, b], [b, b, 0]], periodicity=True)
>>> symm = pychemia.symm.StructureSymmetry(fcc)
>>> symm.number()
225
>>> symm.symbol()
u'Fm-3m'
        """
        assert structure.is_crystal
        assert structure.is_perfect
        self.structure = structure
        # Spglib"s convention for the lattice definition is the transpose of cell
        self._transposed_cell = structure.cell.transpose().copy()
        # Spglib requires numpy floats.
        self._transposed_cell = np.array(self._transposed_cell, dtype='double', order='C')
        self._reduced = np.array(structure.reduced, dtype='double', order='C')

        # Get a list of indices for each atom in structure
        # Indices starting in 1
        self._numbers = np.array([structure.species.index(x) + 1 for x in structure.symbols], dtype='intc')

    @property
    def transposed(self):
        return self._transposed_cell

    @property
    def reduced(self):
        return self._reduced

    @property
    def numbers(self):
        return self._numbers

    def get_symmetry_dataset(self, symprec=1e-5, angle_tolerance=-1.0):

        keys = ('number',
                'hall_number',
                'international',
                'hall',
                'transformation_matrix',
                'origin_shift',
                'rotations',
                'translations',
                'wyckoffs',
                'equivalent_atoms',
                'brv_lattice',
                'brv_types',
                'brv_positions')
        dataset = {}
        for key, data in zip(keys, spg.spglib.spg.dataset(self._transposed_cell,
                                                          self._reduced,
                                                          self._numbers,
                                                          symprec,
                                                          angle_tolerance)):
            dataset[key] = data

        dataset['international'] = dataset['international'].strip()
        dataset['hall'] = dataset['hall'].strip()
        dataset['transformation_matrix'] = np.array(
            dataset['transformation_matrix'], dtype='double', order='C')
        dataset['origin_shift'] = np.array(dataset['origin_shift'], dtype='double')
        dataset['rotations'] = np.array(dataset['rotations'],
                                        dtype='intc', order='C')
        dataset['translations'] = np.array(dataset['translations'],
                                           dtype='double', order='C')
        letters = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ"
        dataset['wyckoffs'] = [letters[x] for x in dataset['wyckoffs']]
        dataset['equivalent_atoms'] = np.array(dataset['equivalent_atoms'],
                                               dtype='intc')
        dataset['brv_lattice'] = np.array(np.transpose(dataset['brv_lattice']),
                                          dtype='double', order='C')
        dataset['brv_types'] = np.array(dataset['brv_types'], dtype='intc')
        dataset['brv_positions'] = np.array(dataset['brv_positions'],
                                            dtype='double', order='C')

        return dataset

    def get_spacegroup(self, symprec=1e-5, angle_tolerance=-1.0):
        """
        Return space group in international table symbol and number
        as a string.

        :param symprec: (float) Symmetry precision
        :param angle_tolerance: (float) angle tolerance for spglib internal routine
        :param symbol_type: (int)
        :return:
        """

        dataset = self.get_symmetry_dataset(symprec=symprec, angle_tolerance=angle_tolerance)
        return "%s (%d)" % (dataset['international'], dataset['number'])

    def spacegroup(self, symprec=1e-5, angle_tolerance=-1.0):
        """
        Computes the space group for the structure with a given
        precision in distances (symprec) and angle tolerance in degrees (angle_tolerance)

        :param symprec: (float) Tolerance of distance between atomic positions and between lengths of lattice vectors
        to be tolerated in the symmetry finding.
        :param angle_tolerance: (float) Tolerance of angle between lattice vectors in degrees to be tolerated in the
        symmetry finding.

        :return: The space group with the symbol and number as a string
        :rtype : (str)
        """

        return self.get_spacegroup(symprec, angle_tolerance)

    def symbol(self, symprec=1e-5, angle_tolerance=-1.0):
        """
        Computes the space group symbol for the structure with a given
        precision in distances (symprec) and angle tolerance in degrees (angle_tolerance)

        :param symprec: (float) Tolerance of distance between atomic positions and between lengths of lattice vectors
        to be tolerated in the symmetry finding.
        :param angle_tolerance: (float) Tolerance of angle between lattice vectors in degrees to be tolerated in the
        symmetry finding.

        :return: The space group symbol as a string
        :rtype : (str)
        """
        return self.spacegroup(symprec, angle_tolerance).split()[0]

    def number(self, symprec=1e-5, angle_tolerance=-1.0):
        """
        Computes the space group number for the structure with a given
        precision in distances (symprec) and angle tolerance in degrees (angle_tolerance)

        :param symprec: (float) Tolerance of distance between atomic positions and between lengths of lattice vectors
        to be tolerated in the symmetry finding.
        :param angle_tolerance: (float) Tolerance of angle between lattice vectors in degrees to be tolerated in the
        symmetry finding.

        :return: The space group number
        :rtype : (int)
        """
        return int(self.spacegroup(symprec, angle_tolerance).split()[-1].strip()[1:-1])

    def refine_cell(self, symprec=1e-5, angle_tolerance=-1.0):
        """
        Refine a pychemia Structure using the tolerances and return a new structure in a Bravais lattice

        :param symprec: (float) Tolerance of distance between atomic positions and between lengths of lattice vectors
                        to be tolerated in the symmetry finding.
        :param angle_tolerance: (float) Tolerance of angle between lattice vectors in degrees to be tolerated in the
                                symmetry finding.
        :return: A new pychemia Structure in a Bravais lattice
        :rtype : (pychemia.Structure)

        Example:

>>> import pychemia
>>> a = 4.05
>>> b = a/2
>>> fcc = pychemia.Structure(symbols=['Au'],
...       cell=[[0, b+1E-5, b-1E-5], [b+1E-5, 0, b-1E-5], [b+1E-5, b-1E-5, 0]], periodicity=True)
>>> symm = pychemia.symm.StructureSymmetry(fcc)
>>> symm.number()
2
>>> symm.symbol()
u'P-1'
>>> fcc2 = symm.refine_cell(symprec=1E-3)
>>> symm2 = pychemia.symm.StructureSymmetry(fcc2)
>>> symm2.number()
225
>>> symm2.symbol()
u'Fm-3m'
        """
        natom = self.structure.natom
        cell = np.array(self.structure.cell.transpose(), dtype='double', order='C')

        pos = np.zeros((natom * 4, 3), dtype='double')
        pos[:natom] = self._reduced

        numbers = np.zeros(natom * 4, dtype='intc')
        numbers[:natom] = np.array(self._numbers, dtype='intc')

        natom_bravais = spg.spglib.spg.refine_cell(cell, pos, numbers, natom, symprec, angle_tolerance)
        if natom_bravais == 0:
            return self.structure.copy()

        cell = np.array(cell.T, dtype='double', order='C')
        reduced = np.array(pos[:natom_bravais], dtype='double', order='C')
        symbols = [self.structure.species[x] for x in (numbers[:natom_bravais] - 1)]

        return Structure(cell=cell, symbols=symbols, reduced=reduced)

    def find_primitive(self, symprec=1e-5, angle_tolerance=1.0):
        """
        Search the primitive pychemia Structure using the tolerances.
        If no primitive cell is found a copy of the original structure is returned

        :param symprec: (float) Tolerance of distance between atomic positions and between lengths of lattice vectors
        to be tolerated in the symmetry finding.
        :param angle_tolerance: (float) Tolerance of angle between lattice vectors in degrees to be tolerated in the
        symmetry finding.

        :return: A new pychemia Structure in a Bravais lattice
        :rtype : (pychemia.Structure)
        """
        # Create copies of the arguments
        cell = np.array(self._transposed_cell, dtype='double', order='C')
        reduced = np.array(self._reduced, dtype='double', order='C')
        numbers = np.array(self._numbers, dtype='intc')

        natom_prim = spg.spglib.spg.primitive(cell, reduced, numbers, symprec, angle_tolerance)
        symbols = [self.structure.species[x] for x in (numbers[:natom_prim] - 1)]
        reduced = reduced[:natom_prim]

        if natom_prim > 0:
            return Structure(cell=cell.T, reduced=reduced, symbols=symbols)
        else:
            return self.structure.copy()

    def crystal_system(self, symprec=1e-5, angle_tolerance=-1.0):

        num = self.number(symprec, angle_tolerance)

        if num < 3:
            return 'Triclinic'
        elif num < 16:
            return 'Monoclinic'
        elif num < 75:
            return 'Orthorhombic'
        elif num < 143:
            return 'Tetragonal'
        elif num < 168:
            return 'Trigonal'
        elif num < 195:
            return 'Hexagonal'
        else:
            return 'Cubic'


def symmetrize(structure, initial_symprec=0.01, final_symprec=0.1, delta_symprec=0.01):
    if structure.natom == 1:
        return structure.copy()
    sym = StructureSymmetry(structure)
    prec = initial_symprec
    while prec < final_symprec:
        if sym.number(symprec=prec) > sym.number(symprec=initial_symprec):
            break
        else:
            prec += delta_symprec
    if prec > final_symprec:
        prec = final_symprec
    new_bravais = sym.refine_cell(symprec=prec)
    sym2 = StructureSymmetry(new_bravais)
    return sym2.find_primitive()
