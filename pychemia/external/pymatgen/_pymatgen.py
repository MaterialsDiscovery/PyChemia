"""
Routines used to interact between pymatgen objects and
pymatdis objects
"""
import pymatgen.io.cifio

import pychemia


def cif2structure(filename, primitive=False, verbose=False):
    """
    Read a given file  as a cif file, convert it into a
    pymatgen structure and finally into a pychemia
    structure
    """
    cifp = pymatgen.io.cifio.CifParser(filename)
    pmg_structs = cifp.get_structures(primitive)
    if len(pmg_structs) != 1:
        if verbose:
            print filename, ' has ', len(pmg_structs), ' structures'
        return None
    ret = pymatgen2pychemia(pmg_structs[0])
    return ret


def pymatgen2pychemia(pmg_struct):
    """
    Converts a pymatgen structure object into a pychemia
    structure object
    """

    #if not all([ x.species_and_occu.is_element for x in pmg_struct]):
    if not pmg_struct.is_ordered:
        #print 'Error: Structure is not Ordered'
        return None

    cell = pmg_struct.lattice_vectors()
    an = pmg_struct.atomic_numbers
    symbols = pychemia.utils.periodic.atomic_symbol(an)
    positions = pmg_struct.cart_coords
    return pychemia.geometry.Structure(cell=cell, positions=positions, symbols=symbols)


def pychemia2pymatgen(structure):
    """
    Converts an pychemia structure into a pymatgen structure object
    """
    lattice = structure.cell
    coords = structure.reduced
    species = structure.symbols
    return pymatgen.Structure(lattice, species, coords)
