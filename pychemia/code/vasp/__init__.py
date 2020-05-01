"""
Module for VASP

The Vienna Ab initio Simulation Package (VASP) is a computer program for atomic scale materials modelling,
e.g. electronic structure calculations and quantum-mechanical molecular dynamics, from first principles.

https://www.vasp.at/
"""

from .kpoints import read_kpoints, write_kpoints
from .poscar import read_poscar, write_poscar, write_potcar, get_potcar_info
from .incar import read_incar, write_incar
from .outcar import VaspOutput, read_vasp_stdout
from .vasp import VaspJob, VaspAnalyser
from .doscar import VaspDoscar
from .queue import write_from_queue
from . import task
from .input import VaspInput
from .output import VaspOutput
from .vaspxml import VaspXML
from .xml_output import parse_vasprun
