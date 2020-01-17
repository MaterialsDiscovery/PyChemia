import os
import sys
from pychemia import HAS_PYMATGEN, pcm_log
from .structure import Structure
from pychemia.code.vasp import read_poscar
from pychemia.code.abinit import AbinitInput


def structure_from_file(structure_file):
    """
    Attempts to reconstruct a PyChemia Structure from the contents of any given file. Valid entries

    :param structure_file: The path to a file where the structure can be reconstructed
    :type structure_file: str
    :return: PyChemia Structure if succeed, None otherwise
    """
    st = None
    basename = os.path.basename(structure_file)
    if not os.path.isfile(structure_file):
        raise ValueError("ERROR: Could not open file '%s'" % structure_file)
    if basename[-4:].lower() == 'json':
        st = Structure.load_json(structure_file)
    elif basename[-3:].lower() == 'cif' and HAS_PYMATGEN:
        import pychemia.external.pymatgen
        st = pychemia.external.pymatgen.cif2structure(structure_file)[0]
    elif 'poscar' in basename.lower():
        st = read_poscar(structure_file)
    elif 'contcar' in basename.lower():
        st = read_poscar(structure_file)
    elif 'abinit' in basename.lower():
        av = AbinitInput(structure_file)
        st = av.get_structure()
    else:
        try:
            st = read_poscar(structure_file)
        except ValueError:
            raise ValueError('Ćould not convert file as POSCAR')
    if st is None:
        pcm_log.debug("ERROR: Could not extract structure from file '%s'" % structure_file)
    return st
