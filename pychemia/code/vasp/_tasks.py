__author__ = 'Guillermo Avendano Franco'

import os
import shutil
import numpy as np
import pychemia
from _incar import InputVariables
from _poscar import read_poscar
from bson.objectid import ObjectId
import logging.handlers


class RelaxPopulation():

    def __init__(self, population, basedir):
        self.population = population
        self.basedir = basedir
        self.VaspJobs = {}

    def create_dirs(self, clean=False):
        if not os.path.isdir(self.basedir):
            os.makedirs(self.basedir)
        elif clean:
            for i in os.listdir(self.basedir):
                shutil.rmtree(self.basedir+os.sep+i)
        for i in self.population.pcdb.entries.find():
            name = self.basedir+os.sep+str(i['_id'])
            if not os.path.isdir(name):
                os.mkdir(name)

    def create_inputs(self):
        kpoints = pychemia.dft.KPoints(kmode='gamma', grid=[4, 4, 4])
        for i in self.population.pcdb.entries.find():
            name = str(i['_id'])
            workdir = self.basedir+os.sep+name
            struct = pychemia.Structure().fromdict(i)
            vj = pychemia.code.vasp.VaspJob(struct, workdir)
            vj.set_kpoints(kpoints)
            inp = pychemia.code.vasp.InputVariables()
            inp.set_minimum()
            #inp.set_break_conditions()
            inp.set_ion_relax()
            inp.set_rough_relaxation()
            vj.set_input_variables(inp)
            vj.write_all()
            vj.input_variables.set_encut(ENCUT=1.1, POTCAR=workdir+os.sep+'POTCAR')
            vj.write_incar()
            self.VaspJobs[name] = vj

    def update(self, workdir):

        working_id = os.path.basename(workdir)
        vj = self.VaspJobs[working_id]
        mongoid = ObjectId(working_id)

        if not os.path.isfile(workdir+os.sep+'PROCAR'):
            return True
        else:
            if os.path.isfile(workdir+os.sep+'OUTCAR'):
                print '-'
                vo = pychemia.code.vasp.VaspOutput(workdir+os.sep+'OUTCAR')
                info = vo.relaxation_info()
                if len(info) != 3:
                    print '['+str(working_id)+']'+' Missing some data in OUTCAR (forces or stress)'
                    return True

                for i in info:
                    print '['+str(working_id)+'] %20s %12.5e' % (i, info[i])

                if info['avg_force'] < 0.001:
                    if info['avg_stress_diag'] < 0.001:
                        if info['avg_stress_non_diag'] < 0.001:
                            wf = open(workdir+os.sep+'COMPLETE', 'w')
                            for i in info:
                                wf.write("%15s %12.3f" % (i, info[i]))
                            wf.close()
                            return False

                if info['avg_force'] < 0.01:
                    if info['avg_stress_diag'] < 0.01:
                        if info['avg_stress_non_diag'] < 0.01:
                            vj.input_variables.variables['ISIF'] = 3
                        else:
                            vj.input_variables.variables['ISIF'] = 5
                    else:
                        vj.input_variables.variables['ISIF'] = 7
                else:
                    vj.input_variables.variables['ISIF'] = 2
                print '['+str(working_id)+']'+'ISIF   : ', vj.input_variables.variables['ISIF']

                if vj.input_variables.variables['EDIFF'] > 1E-7:
                    vj.input_variables.variables['EDIFF'] /= 10
                print '['+str(working_id)+']'+'EDIFF  : ', vj.input_variables.variables['EDIFF']

                if vj.input_variables.variables['EDIFFG'] < -1E-5:
                    vj.input_variables.variables['EDIFFG'] /= 10
                print '['+str(working_id)+']'+'EDIFFG : ', vj.input_variables.variables['EDIFFG']
                print '-'

                for i in ['INCAR', 'POSCAR', 'OUTCAR', 'vasp.stdout', 'vasp.stderr']:
                    if os.path.exists(workdir+os.sep+i):
                        log = logging.handlers.RotatingFileHandler(workdir+os.sep+i, maxBytes=1, backupCount=1000)
                        log.doRollOver()

                vj.structure = read_poscar(workdir+os.sep+'CONTCAR')
                vj.write_all()

                self.population.update_entry()
                vj.save_json(workdir+os.sep+'PyChemia.entry')

                return True
            else:
                return True

    @property
    def workdirs(self):
        return [self.basedir+os.sep+name for name in self.population.entries_ids]


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
        return success

    def postprocess(self):
        pass
