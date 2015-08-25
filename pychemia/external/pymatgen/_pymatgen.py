"""
Routines used to interact between pymatgen objects and
pymatdis objects
"""
import numpy as np

try:
    import pymatgen.io.cif
except ImportError:
    raise ImportError
from pychemia import Structure
from pymatgen import Structure as PMG_Structure


def cif2structure(filename, primitive=False):
    """
    Read a given file  as a cif file, convert it into a
    pymatgen structure and finally into a pychemia
    structure
    """
    cifp = pymatgen.io.cif.CifParser(filename)
    pmg_structs = cifp.get_structures(primitive)
    ret = []
    for i in pmg_structs:
        ret.append(pymatgen2pychemia(i))
    return ret


def pymatgen2pychemia(pmg_struct):
    """
    Converts a pymatgen structure object into a pychemia
    structure object
    """

    # if not pmg_struct.is_ordered:
    #     return None

    cell = pmg_struct.lattice_vectors()
    reduced = np.array([])
    sites = []
    symbols = []
    occupancies = []
    index = 0
    for i in pmg_struct:
        for j in i.species_and_occu:
            symbols.append(j.symbol)
            occupancies.append(i.species_and_occu[j])
            sites.append(index)
        reduced = np.concatenate((reduced, i.frac_coords))
        index += 1

    reduced.reshape((-1, 3))
    return Structure(cell=cell, reduced=reduced, symbols=symbols, sites=sites, occupancies=occupancies)


def pychemia2pymatgen(structure):
    """
    Converts an pychemia structure into a pymatgen structure object
    """
    lattice = structure.cell
    coords = structure.reduced
    species = structure.symbols
    return PMG_Structure(lattice, species, coords)
