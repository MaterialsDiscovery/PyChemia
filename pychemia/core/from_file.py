__author__ = 'Guillermo Avendano-Franco'

import os
import sys
from .structure import Structure

def structure_from_file(structure_file):

    from pychemia.code.vasp import read_poscar
    from pychemia.code.abinit import InputVariables

    st = None
    if not os.path.isfile(structure_file):
        print "ERROR: Could not open file '%s'" % structure_file
        sys.exit(1)
    if structure_file[-4:].lower() == 'json':
        st = Structure.load_json(structure_file)
    elif structure_file[-3:].lower() == 'cif' and pychemia.HAS_PYMATGEN:
        import pychemia.external.pymatgen
        st = pychemia.external.pymatgen.cif2structure(structure_file)[0]
    elif structure_file[-6:].lower() == 'poscar':
        st = read_poscar(structure_file)
    elif structure_file[-6:].lower() == 'contcar':
        st = read_poscar(structure_file)
    elif structure_file[-6:].lower() == 'abinit.in':
        av = InputVariables(structure_file)
        av.get_structure()
    else:
        try:
            st = read_poscar(structure_file)
        except ValueError:
            print "ERROR: Could not extract structure from file '%s'" % structure_file
            exit(1)
    return st
