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

from pyspglib import spglib

import ase.atoms
import ase.io

import pychemia


def cif2structure(filename, primitive=False, symprec=0.001):
    """
    Read a given file  as a cif file, convert it into a
    ase atoms object and finally into a pychemia
    structure
    """
    aseatoms = ase.io.read(filename)
    if primitive:
        #lattice, scaled_positions, numbers = spglib.refine_cell(aseatoms, symprec)
        #ref_atoms=ase.atoms(cell=lattice,scaled_positions=scaled_positions,numbers=numbers)
        lattice, scaled_positions, numbers = spglib.find_primitive(aseatoms, symprec)
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
    return pychemia.core.Structure(cell=cell, positions=positions, symbols=symbols)


def pychemia2ase(structure):
    """
    Converts an pychemia structure into a ase atoms object
    """
    cell = structure.cell
    scaled_positions = structure.reduced
    symbols = structure.symbols
    return ase.atoms.Atoms(cell=cell, scaled_positions=scaled_positions, symbols=symbols)
