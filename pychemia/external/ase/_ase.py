import pychemia

try:
    import spglib as spg
    USE_SPGLIB = True
except ImportError:
    USE_SPGLIB = False
    # from pyspglib import spglib as spg

"""
Interface with the Atomic Simulation Environment (ASE)
Convert ase Atoms objects into pymatdis crystal objects
"""

__author__ = "Guillermo Avendano-Franco"
__copyright__ = "Copyright 2014"
__version__ = "0.1"
__maintainer__ = "Guillermo Avendano-Franco"
__email__ = "gtux.gaf@gmail.com"
__status__ = "Development"
__date__ = "March 31, 2014"

try:
    import ase.atoms
    import ase.io
except ImportError:
    print('The module ase could not be imported')
    raise ImportError


def cif2structure(filename, primitive=False, symprec=0.001):
    """
    Read a given file  as a cif file, convert it into a
    ase atoms object and finally into a pychemia
    structure

    :param filename: (str) Filename of the CIF file that will be read
    :param primitive: (boolean) if primitive cell should be returned
    :param symprec: (float) Desired precision for computing the primitive cell
    :return:
    """
    aseatoms = ase.io.read(filename)
    if primitive and USE_SPGLIB:
        lattice, scaled_positions, numbers = spg.find_primitive(aseatoms, symprec)
        if lattice is not None and scaled_positions is not None and numbers is not None:
            fin_atoms = ase.atoms.Atoms(cell=lattice, scaled_positions=scaled_positions, numbers=numbers, pbc=True)
        else:
            fin_atoms = aseatoms
    else:
        fin_atoms = aseatoms
    ret = ase2pychemia(fin_atoms)
    return ret


def ase2pychemia(aseatoms):
    """
    Converts an ase atoms object into a pychemia
    structure object
    """
    cell = aseatoms.get_cell()
    an = aseatoms.get_atomic_numbers()
    symbols = pychemia.utils.periodic.atomic_symbol(an)
    positions = aseatoms.get_positions()
    return pychemia.Structure(cell=cell, positions=positions, symbols=symbols)


def pychemia2ase(structure):
    """
    Converts an pychemia structure into a ase atoms object

    :param structure: (pychemia.Structure) PyChemia Structure to convert into a ASE Atoms object
    :return:
    """
    cell = structure.cell
    scaled_positions = structure.reduced
    symbols = structure.symbols
    return ase.atoms.Atoms(cell=cell, scaled_positions=scaled_positions, symbols=symbols)
