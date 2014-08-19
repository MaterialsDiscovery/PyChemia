__author__ = 'Guillermo Avendano-Franco'

import numpy as np

try:
    import pyspglib._spglib as spg
except ImportError:
    msg = """Spglib required. Please install pyspglib from spglib.
    Try to test that you are able to import correctly the module
    using:

    import pyspglib._spglib as spg
    """
    raise ImportError(msg)


class Symmetry(object):
    """
    Takes a pychemia.core.Structure object and compute the space group
    of that structure.
    Uses pyspglib to perform various symmetry finding operations.
    Two parameters affect the algorithm, the default values are tuned for
    structures optimized using ab-initio codes

    Args:
        structure (pychemia.core.Structure): Structure to find symmetry
        tolerance (float): Defaults to 1e-1,
        angle_tolerance (float): Angle tolerance. Defaults to 5
    """

    def __init__(self, structure, tolerance=1e-1, angle_tolerance=5):
        self._tolerance = tolerance
        self._angle_tolerance = angle_tolerance
        self._structure = structure
        #Spglib"s convention for the lattice definition is the transpose of cell
        self._transposed_cell = structure.get_cell().cell.transpose()
        #Spglib requires numpy floats.
        self._transposed_cell = np.array(
            self._transposed_cell, dtype='double', order='C')
        self._reduced = np.array(structure.reduced, dtype='double', order='C')

        numbers = list(structure.symbols)
        index = 0
        for i in range(structure.natom):
            index += 1
            if numbers[i].isalpha():
                for j in range(i, structure.natom):
                    if numbers[j] == structure.symbols[i]:
                        numbers[j] = str(index)
            else:
                index -= 1

        for i in range(structure.natom):
            numbers[i] = int(numbers[i])

        self._numbers = np.array(numbers, dtype='intc')
        self._spacegroup_data = spg.spacegroup(
            self._transposed_cell.copy(), self._reduced.copy(),
            self._numbers, self._tolerance, self._angle_tolerance)

    @property
    def symbol(self):
        """
        Get the spacegroup symbol string.

        Returns:
            (str): Spacegroup symbol as string.
        """
        return self._spacegroup_data.split()[0]

    @property
    def number(self):
        """
        Get the international spacegroup number.

        Returns:
            (int): International spacegroup number.
        """
        return int(self._spacegroup_data.split()[-1].strip()[1:-1])
