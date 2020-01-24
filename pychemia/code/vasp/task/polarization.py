from ...tasks import Task
import os
import shutil
import numpy as np
from ..kpoints import write_kpoints
from ..poscar import write_poscar
from ..incar import write_incar
from ..input import VaspInput
from ....runner import Runner

__author__ = 'Guillermo Avendano-Franco'


class Polarization(Task):
    def report(self):
        pass

    def save(self, filename=None):
        pass

    def plot(self):
        pass

    def status(self):
        pass

    def load(self):
        pass

    def __init__(self, structure, workdir, potcar_filepath, external=None, maxfield=0, stepfield=0, executable='vasp'):

        self.external = None
        self.potcar = potcar_filepath

        if external is not None and external in ['electric', 'magnetic']:
            self.external = external

        self.maxfield = abs(maxfield)
        self.stepfield = abs(stepfield)
        task_params = {'maxfield': self.maxfield, 'stepfield': self.stepfield, 'external': self.external,
                       'portcar': self.potcar}
        Task.__init__(self, structure=structure, task_params=task_params, workdir=workdir, executable=executable)

    def initialize(self, kpoints, cleandir=False):
        if not os.path.isdir(self.workdir):
            os.mkdir(self.workdir)
        elif cleandir:
            shutil.rmtree(self.workdir)
            os.mkdir(self.workdir)

        if self.maxfield > 0:
            signs = ['+', '-']
            fields = np.arange(0, self.maxfield, self.stepfield)
        else:
            signs = ['']
            fields = [0]

        for sign in signs:
            for field in fields:
                pathname = self.workdir + os.sep + 'Field' + sign + str(field)
                if not os.path.isdir(pathname):
                    os.mkdir(pathname)
                for step in ['SCF', 'BANDS', 'berry_IGPAR1', 'berry_IGPAR2', 'berry_IGPAR3']:
                    if not os.path.isdir(pathname + os.sep + step):
                        os.mkdir(pathname + os.sep + step)
                    write_kpoints(kpoints, pathname + os.sep + step + os.sep + 'KPOINTS')
                    if not os.path.lexists(pathname + os.sep + step + os.sep + 'POTCAR'):
                        os.symlink(os.path.abspath(self.potcar), pathname + os.sep + step + os.sep + 'POTCAR')

    def run(self, mode='local'):

        if self.maxfield > 0:
            signs = ['+', '-']
            fields = np.arange(0, self.maxfield, self.stepfield)
        else:
            signs = ['']
            fields = [0]

        success = True
        for sign in signs:
            for field in fields:
                pathname = self.workdir + os.sep + 'Field' + sign + str(field)
                for step in ['SCF', 'BANDS', 'berry_IGPAR1', 'berry_IGPAR2', 'berry_IGPAR3']:
                    tk = Task()
                    tk.minimum(ENCUT=200, POTCAR=self.potcar)
                    tk.vaspinput['EDIFF'] = 1E-9
                    tk.vaspinput['IBRION'] = -1
                    tk.vaspinput['NSW'] = 0
                    if step == 'SCF':
                        iv = VaspInput(variables=tk.vaspinput)
                        write_incar(iv, pathname + os.sep + step + os.sep + 'INCAR')
                        write_poscar(structure=self.structure, filepath=pathname + os.sep + step + os.sep + 'POSCAR')
                        options = {'nproc': 4, 'code_bin': '/home/guilleaf/local/src/vasp.5.3/vasp', 'mpi': True}
                        runner = Runner('vasp', 'local', options)
                        runner.run(dirpath=pathname + os.sep + step)
                        rf = open(pathname + os.sep + step + os.sep + 'vasp.stdout')
                        if rf.readlines()[-2].strip() != 'writing wavefunctions':
                            print('Finished writing the Wavefunctions')
                            success = True
                        else:
                            success = False
                    if step == 'BANDS':
                        shutil.copy('KPOINTS-path', pathname + os.sep + step + os.sep + 'KPOINTS')
                        shutil.copy2(pathname + os.sep + 'SCF/CHGCAR', pathname + os.sep + step)
                        tk.vaspinput['ICHARG'] = 11
                        iv = VaspInput(variables=tk.vaspinput)
                        write_incar(iv, pathname + os.sep + step + os.sep + 'INCAR')
                        write_poscar(structure=self.structure, filepath=pathname + os.sep + step + os.sep + 'POSCAR')
                        options = {'nproc': 4, 'code_bin': 'vasp', 'mpi': True}
                        runner = Runner('vasp', 'local', options)
                        runner.run(dirpath=pathname + os.sep + step)
                        rf = open(pathname + os.sep + step + os.sep + 'vasp.stdout')
                        if rf.readlines()[-2].strip() != 'writing wavefunctions':
                            print('Finished writing the Wavefunctions')
                            success = True
                        else:
                            success = False
                    if step.startswith('berry_IGPAR'):
                        shutil.copy2(pathname + os.sep + 'SCF/CHGCAR', pathname + os.sep + step)
                        tk.vaspinput['ICHARG'] = 1
                        tk.vaspinput['EDIFF'] = 1E-9
                        tk.vaspinput['IBRION'] = -1
                        tk.vaspinput['NSW'] = 0
                        tk.vaspinput['LBERRY'] = True
                        tk.vaspinput['NPPSTR'] = 12
                        tk.vaspinput['IGPAR'] = int(step[-1])
                        if 'ISIF' in tk.vaspinput:
                            tk.vaspinput.pop('ISIF')

                        iv = VaspInput(variables=tk.vaspinput)
                        write_incar(iv, pathname + os.sep + step + os.sep + 'INCAR')
                        write_poscar(structure=self.structure, filepath=pathname + os.sep + step + os.sep + 'POSCAR')
                        options = {'nproc': 4, 'code_bin': '/home/guilleaf/local/src/vasp.5.3/vasp', 'mpi': True}
                        runner = Runner('vasp', 'local', options)
                        runner.run(dirpath=pathname + os.sep + step)
                        rf = open(pathname + os.sep + step + os.sep + 'vasp.stdout')
                        if rf.readlines()[-2].strip() == 'writing wavefunctions':
                            print('Finished writing the Wavefunctions')
                            success = True
                        else:
                            success = False
        return success

    def postprocess(self):
        pass
