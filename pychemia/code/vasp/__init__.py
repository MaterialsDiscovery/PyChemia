"""
Module for VASP

The Vienna Ab initio Simulation Package (VASP) is a computer program for atomic scale materials modelling,
e.g. electronic structure calculations and quantum-mechanical molecular dynamics, from first principles.

https://www.vasp.at/
"""

from ._kpoints import read_kpoints, write_kpoints
from ._poscar import read_poscar, write_poscar, write_potcar
from ._incar import read_incar, write_incar, InputVariables
from ._outcar import VaspOutput, read_vasp_stdout
from ._vasp import VaspJob, VaspAnalyser
from ._doscar import VaspDoscar
from ._queue import write_from_queue
from . import task
