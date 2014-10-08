__author__ = 'Guillermo Avendano Franco'

import os
import shutil
import numpy as np
from pychemia.code.vasp._kpoints import save_KPOINTS
from pychemia.code.vasp._incar import save_INCAR, InputVariables
from pychemia.code.vasp._poscar import save_POSCAR
import pychemia


class Tasks():

    def __init__(self, vaspinput=None):
        """
        Create a VASP tasks object completing information in the input for running
        specific calculations using VASP

        :param input:
        """
        if vaspinput is not None:
            self.vaspinput = vaspinput
        else:
            self.vaspinput = {}

    def minimum(self, PREC='Normal', ISPIN=2, ENCUT=300, LREAL=False, ISMEAR=0, LORBIT=11, POTCAR=None):
        self.vaspinput['PREC'] = PREC
        self.vaspinput['ENCUT'] = ENCUT
        self.vaspinput['LREAL'] = LREAL
        self.vaspinput['ISMEAR'] = ISMEAR
        self.vaspinput['ISPIN'] = ISPIN
        self.vaspinput['LORBIT'] = LORBIT
        if POTCAR is not None and ENCUT < 10:
            maxvalue = 0
            if not os.path.isfile(POTCAR):
                raise ValueError('Not such file', POTCAR)
            rf = open(POTCAR)
            for line in rf.readlines():
                if 'ENMAX' in line:
                    list4line = line.split()
                    assert(list4line[0].strip() == 'ENMAX')
                    value = list4line[2].strip()
                    if value[-1] == ';':
                        value = value[:-1]
                    value = float(value)
                    if value > maxvalue:
                        maxvalue = value
            rf.close()
            if ENCUT < 10:
                self.vaspinput['ENCUT'] = ENCUT*maxvalue
            else:
                self.vaspinput['ENCUT'] = maxvalue

    def ion_relax(self, NSW=100, ISIF=2):
        self.vaspinput['IBRION'] = 2
        self.vaspinput['NSW'] = NSW
        self.vaspinput['ISIF'] = ISIF

    def ion_tolerances(self, EDIFF='1E-7', EDIFFG='-1E-3'):
        self.vaspinput['EDIFF'] = EDIFF
        self.vaspinput['EDIFFG'] = EDIFFG


class Polarization():

    def __init__(self, structure, path, potcar_filepath, external=None, maxfield=0, stepfield=0):

        self.structure = structure
        self.path = path
        self.external = None
        self.potcar = potcar_filepath

        if external is not None and extenal in ['electric', 'magnetic']:
            self.external = external

        self.maxfield = abs(maxfield)
        self.stepfield = abs(stepfield)

    def initialize(self, kpoints,  cleandir=False):
        if not os.path.isdir(self.path):
            os.mkdir(self.path)
        elif cleandir:
            shutil.rmtree(self.path)
            os.mkdir(self.path)

        if self.maxfield > 0:
            signs = ['+', '-']
            fields = np.arange(0, self.maxfield, self.stepfield)
        else:
            signs = ['']
            fields = [0]

        for sign in signs:
            for field in fields:
                pathname = self.path+os.sep+'Field'+sign+str(field)
                if not os.path.isdir(pathname):
                    os.mkdir(pathname)
                for step in ['SCF', 'BANDS', 'berry_IGPAR1', 'berry_IGPAR2', 'berry_IGPAR3']:
                    if not os.path.isdir(pathname + os.sep + step):
                        os.mkdir(pathname + os.sep + step)
                    save_KPOINTS(kpoints, pathname + os.sep + step+os.sep+'KPOINTS')
                    if not os.path.lexists(pathname + os.sep + step+os.sep+'POTCAR'):
                        os.symlink(os.path.abspath(self.potcar), pathname + os.sep + step+os.sep+'POTCAR')

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
                pathname = self.path+os.sep+'Field'+sign+str(field)
                for step in ['SCF', 'BANDS', 'berry_IGPAR1', 'berry_IGPAR2', 'berry_IGPAR3']:
                    tk = Tasks()
                    tk.minimum(ENCUT=200, POTCAR=self.potcar)
                    tk.vaspinput['EDIFF'] = 1E-9
                    tk.vaspinput['IBRION'] = -1
                    tk.vaspinput['NSW'] = 0
                    if step == 'SCF':
                        iv = InputVariables(variables=tk.vaspinput)
                        save_INCAR(iv, pathname + os.sep + step+os.sep+'INCAR')
                        save_POSCAR(structure=self.structure, filepath=pathname + os.sep + step+os.sep+'POSCAR')
                        options = {'nproc': 4, 'code_bin': '/home/guilleaf/local/src/vasp.5.3/vasp', 'mpi': True}
                        runner = pychemia.runner.Runner('vasp', 'local', options)
                        p = runner.run(dirpath=pathname + os.sep + step)
                        rf = open(pathname + os.sep + step + os.sep + 'vasp.stdout')
                        if rf.readlines()[-2].strip() != 'writing wavefunctions':
                            print 'Finished writing the Wavefunctions'
                            success = True
                        else:
                            success = False
                    if step == 'BANDS':
                        shutil.copy('KPOINTS-path', pathname + os.sep + step + os.sep + 'KPOINTS')
                        shutil.copy2(pathname + os.sep + 'SCF/CHGCAR', pathname + os.sep + step)
                        tk.vaspinput['ICHARG'] = 11
                        iv = InputVariables(variables=tk.vaspinput)
                        save_INCAR(iv, pathname + os.sep + step+os.sep+'INCAR')
                        save_POSCAR(structure=self.structure, filepath=pathname + os.sep + step+os.sep+'POSCAR')
                        options = {'nproc': 4, 'code_bin': '/home/guilleaf/local/src/vasp.5.3/vasp', 'mpi': True}
                        runner = pychemia.runner.Runner('vasp', 'local', options)
                        p = runner.run(dirpath=pathname + os.sep + step)
                        rf = open(pathname + os.sep + step + os.sep + 'vasp.stdout')
                        if rf.readlines()[-2].strip() != 'writing wavefunctions':
                            print 'Finished writing the Wavefunctions'
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

                        iv = InputVariables(variables=tk.vaspinput)
                        save_INCAR(iv, pathname + os.sep + step+os.sep+'INCAR')
                        save_POSCAR(structure=self.structure, filepath=pathname + os.sep + step+os.sep+'POSCAR')
                        options = {'nproc': 4, 'code_bin': '/home/guilleaf/local/src/vasp.5.3/vasp', 'mpi': True}
                        runner = pychemia.runner.Runner('vasp', 'local', options)
                        p = runner.run(dirpath=pathname + os.sep + step)
                        rf = open(pathname + os.sep + step + os.sep + 'vasp.stdout')
                        if rf.readlines()[-2].strip() == 'writing wavefunctions':
                            print 'Finished writing the Wavefunctions'
                            success = True
                        else:
                            success = False

    def postprocess(self):
        pass
