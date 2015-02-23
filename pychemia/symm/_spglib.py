__author__ = 'Guillermo Avendano-Franco'

import numpy as np
from pychemia import Structure

try:
    import pyspglib._spglib as spg
except ImportError:
    msg = """Spglib required. Please install pyspglib from spglib.
    Try to test that you are able to import correctly the module
    using:

    import pyspglib._spglib as spg
    """
    raise ImportError(msg)


class StructureSymmetry(object):
    """
    Takes a pychemia.Structure object and creates an object with methods
    to identify symmetry groups and other symmetry related operations.
    Uses pyspglib to perform various symmetry finding operations.
    """

    def __init__(self, structure):
        """
        Creates a new StructureSymmetry object for a given structure

        :param structure:
        :rtype : (Symmetry)
        """
        assert structure.is_crystal
        assert structure.is_perfect
        self._structure = structure
        # Spglib"s convention for the lattice definition is the transpose of cell
        self._transposed_cell = structure.cell.transpose().copy()
        # Spglib requires numpy floats.
        self._transposed_cell = np.array(self._transposed_cell, dtype='double', order='C')
        self._reduced = np.array(structure.reduced, dtype='double', order='C')

        # Get a list of indices for each atom in structure
        # Indices starting in 1
        self._numbers = np.array([structure.species.index(x)+1 for x in structure.symbols], dtype='intc')

    def spacegroup(self, symprec=1e-1, angle_tolerance=5):
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
        return spg.spacegroup(
            self._transposed_cell.copy(), self._reduced.copy(),
            self._numbers, symprec, angle_tolerance)

    def symbol(self, symprec=1e-1, angle_tolerance=5):
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

    def number(self, symprec=1e-1, angle_tolerance=5):
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

    def refine_cell(self, symprec=1e-1, angle_tolerance=5):
        """
        Refine a pychemia Structure using the tolerances and return a new structure in a Bravais lattice

        :param symprec: (float) Tolerance of distance between atomic positions and between lengths of lattice vectors
        to be tolerated in the symmetry finding.
        :param angle_tolerance: (float) Tolerance of angle between lattice vectors in degrees to be tolerated in the
        symmetry finding.

        :return: A new pychemia Structure in a Bravais lattice
        :rtype : (pychemia.Structure)
        """
        natom = self._structure.natom
        cell = np.array(self._structure.cell.transpose(), dtype='double', order='C')

        pos = np.zeros((natom * 4, 3), dtype='double')
        pos[:natom] = self._reduced

        numbers = np.zeros(natom * 4, dtype='intc')
        numbers[:natom] = np.array(self._numbers, dtype='intc')

        natom_bravais = spg.refine_cell(cell, pos, numbers, natom, symprec, angle_tolerance)
        cell = np.array(cell.T, dtype='double', order='C')
        reduced = np.array(pos[:natom_bravais], dtype='double', order='C')

        symbols = [self._structure.species[x] for x in (numbers[:natom_bravais]-1)]

        return Structure(cell=cell, symbols=symbols, reduced=reduced)

    def find_primitive(self, symprec=1e-1, angle_tolerance=5):
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

        natom_prim = spg.primitive(cell, reduced, numbers, symprec, angle_tolerance)
        symbols = [self._structure.species[x] for x in (numbers[:natom_prim]-1)]
        reduced = reduced[:natom_prim]

        if natom_prim > 0:
            return Structure(cell=cell.T, reduced=reduced, symbols=symbols)
        else:
            return self._structure.copy()


def symmetrize(structure, initial_symprec=0.01, final_symprec=0.2, delta_symprec=0.01):

    sym = StructureSymmetry(structure)
    prec = initial_symprec
    while prec < final_symprec:
        if sym.number(symprec=prec) > sym.number(symprec=initial_symprec):
            break
        else:
            prec += delta_symprec
    new_bravais = sym.refine_cell(symprec=prec)
    sym2 = StructureSymmetry(new_bravais)
    return sym2.find_primitive()
