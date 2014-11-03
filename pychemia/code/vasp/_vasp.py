__author__ = 'guilleaf'

import os
from _poscar import write_poscar, write_potcar, read_poscar
from _kpoints import write_kpoints, read_kpoints
from _incar import write_incar, InputVariables, read_incar
from pychemia import Structure
from pychemia.dft import KPoints
import json


class VaspJob():

    def __init__(self, structure=None, workdir=None):

        self.structure = structure
        self.workdir = workdir
        self.input_variables = None
        self.potcar_setup = None
        self.potcar_pspfiles = None
        self.potcar_pspdir = 'potpaw_PBE'
        self.kpoints = None
        self.outcar = None
        self.poscar_setup = None

    def _check_workdir(self):

        if self.workdir is None:
            raise ValueError("A proper working directory has not been setup")
        elif not os.path.isdir(self.workdir):
            os.makedirs(self.workdir)

    def write_poscar(self):

        self._check_workdir()
        assert(isinstance(self.structure, Structure))

        write_poscar(self.structure, filepath=self.workdir+os.sep+'POSCAR')

    def write_potcar(self):

        self._check_workdir()
        assert(isinstance(self.structure, Structure))

        pspfiles = write_potcar(self.structure, filepath=self.workdir+os.sep+'POTCAR', pspdir=self.potcar_pspdir,
                                options=self.potcar_setup, pspfiles=self.potcar_pspfiles)

        self.potcar_pspfiles = pspfiles

    def set_kpoints(self, kpoints):

        assert isinstance(kpoints, KPoints)
        self.kpoints = kpoints

    def write_kpoints(self):

        self._check_workdir()
        assert(isinstance(self.structure, Structure))

        write_kpoints(self.kpoints, filepath=self.workdir+os.sep+'KPOINTS')

    def write_incar(self):

        self._check_workdir()
        assert(isinstance(self.input_variables, InputVariables))

        write_incar(self.input_variables, filepath=self.workdir+os.sep+'INCAR')

    def set_input_variables(self, input_variables):

        self._check_workdir()
        assert(isinstance(input_variables, InputVariables))
        self.input_variables = input_variables

    def write_all(self):

        self.write_incar()
        self.write_kpoints()
        self.write_poscar()
        self.write_potcar()

    @property
    def variables(self):
        return self.input_variables.variables

    def todict(self):
        ret = {'structure': self.structure.todict(),
               'potcar_pspfiles': self.potcar_pspfiles,
               'potcar_setup': self.potcar_setup,
               'potcar_workdir': self.workdir,
               'variables': self.input_variables.variables,
               'kpoints': self.kpoints.todict(),
               'outcar': self.outcar,
               'poscar_setup': self.poscar_setup}
        return ret

    def save_json(self, filename):

        filep = open(filename, 'w')
        json.dump(self.todict(), filep, sort_keys=True, indent=4, separators=(',', ': '))
        filep.close()

    def read_incar(self):
        self.input_variables = read_incar(self.workdir+os.sep+'INCAR')

    def read_kpoint(self):
        self.kpoints = read_kpoints(self.workdir+os.sep+'KPOINTS')

    def read_poscar(self):
        self.structure = read_poscar(self.workdir+os.sep+'POSCAR')

