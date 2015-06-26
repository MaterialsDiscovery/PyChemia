__author__ = 'Guillermo Avendano Franco'

import os
import time
import shutil
import logging.handlers
import json
import numpy as np

from _incar import InputVariables
from _outcar import VaspOutput, read_vasp_stdout
from _vasp import VaspJob, VaspAnalyser
from _poscar import read_poscar
from _kpoints import read_kpoints
from pychemia.utils.mathematics import round_small
from pychemia.dft import KPoints
from pychemia import pcm_log
from pychemia.code import Relaxator
from pychemia.runner import Runner


class RelaxPopulation:
    def __init__(self, population, basedir, target_force=1E-2, target_stress=1E-2):
        self.population = population
        self.basedir = basedir
        self.vasp_jobs = {}
        self.runs = {}
        self.runner = None
        self.status = {}
        self.target_force = target_force
        self.target_stress = target_stress

    def create_dirs(self, clean=False):
        if not os.path.isdir(self.basedir):
            os.makedirs(self.basedir)
        elif clean:
            for i in os.listdir(self.basedir):
                shutil.rmtree(self.basedir + os.sep + i)
        for i in self.population.pcdb.entries.find():
            name = self.basedir + os.sep + str(i['_id'])
            if not os.path.isdir(name):
                os.mkdir(name)

    def create_inputs(self, density_of_kpoints=10000, encut=1.0):
        # kpoints = KPoints(kmode='gamma', grid=[4, 4, 4])
        kpoints = KPoints()
        for entry in self.population.pcdb.entries.find():
            name = str(entry['_id'])
            workdir = self.basedir + os.sep + name
            structure = self.population.db.get_structure(entry['_id'])
            kpoints.set_optimized_grid(structure.lattice, density_of_kpoints=density_of_kpoints)
            print kpoints
            vj = VaspJob()
            vj.initialize(workdir=workdir, structure=structure, kpoints=kpoints)
            inp = InputVariables()
            inp.set_rough_relaxation()
            vj.set_input_variables(inp)
            vj.write_potcar()
            vj.input_variables.set_encut(ENCUT=encut, POTCAR=workdir + os.sep + 'POTCAR')
            vj.set_inputs()
            self.vasp_jobs[name] = vj
            self.runs[name] = 0
            self.status[name] = ['ACTIVE']

    def add_status(self, entry_id, value):
        if value not in self.status[entry_id]:
            self.status[entry_id].append(value)

    def del_status(self, entry_id, value):
        if value in self.status[entry_id]:
            self.status[entry_id].remove(value)

    def flip_status(self, entry_id, oldvalue, newvalue):
        self.del_status(entry_id, oldvalue)
        self.add_status(entry_id, newvalue)

    def modify_input(self, entry_id):
        if 'RELAXED' not in self.status[entry_id] and 'NOPROCAR' not in self.status[entry_id] \
                and 'NOOUTCAR' not in self.status[entry_id]:
            return True
        else:
            return False

    def update(self, workdir):
        """
        This routine determines how to proceed with the relaxation
        for one specific work directory

        :param workdir: (str) String representation of the id in the mongodb
        :return:
        """

        # workdir = self.basedir + os.sep + entry_id
        entry_id = os.path.basename(workdir)
        vj = self.vasp_jobs[entry_id]
        runj = self.runs[entry_id]
        if os.path.isfile(workdir + os.sep + 'OUTCAR'):
            vj.get_outputs()
        self.update_history(entry_id)

        if os.path.isfile(workdir + os.sep + 'RELAXED'):
            self.add_status(entry_id, 'RELAXED')
        elif not os.path.isfile(workdir + os.sep + 'PROCAR'):
            self.add_status(entry_id, 'NOPROCAR')
        else:
            self.del_status(entry_id, 'NOPROCAR')
            if not os.path.isfile(workdir + os.sep + 'OUTCAR'):
                self.add_status(entry_id, 'NOOUTCAR')
            else:
                self.del_status(entry_id, 'NOOUTCAR')
                print '-'
                vo = VaspOutput(workdir + os.sep + 'OUTCAR')
                relaxation_info = vo.relaxation_info()
                if len(relaxation_info) != 3:
                    print '[' + str(entry_id) + ']' + ' Missing some data in OUTCAR (forces or stress)'
                    self.add_status(entry_id, 'NOOUTCAR')

                print '[' + str(entry_id) + ']' + 'Results:'
                for i in relaxation_info:
                    print '[' + str(entry_id) + '] %20s %12.5e' % (i, relaxation_info[i])

                # Conditions to consider the structure relaxed
                if relaxation_info['avg_force'] < self.target_force:
                    if relaxation_info['avg_stress_diag'] < self.target_stress:
                        if relaxation_info['avg_stress_non_diag'] < self.target_stress:
                            wf = open(workdir + os.sep + 'RELAXED', 'w')
                            for i in relaxation_info:
                                wf.write("%15s %12.3f" % (i, relaxation_info[i]))
                            wf.close()
                            wf = open(workdir + os.sep + 'COMPLETE', 'w')
                            for i in relaxation_info:
                                wf.write("%15s %12.3f" % (i, relaxation_info[i]))
                            wf.close()
                            self.add_status(entry_id, 'RELAXED')

        if self.modify_input(entry_id):

            # How to change ISIF
            if relaxation_info['avg_force'] < 0.1:
                if relaxation_info['avg_stress_diag'] < 0.1:
                    if relaxation_info['avg_stress_non_diag'] < 0.1:
                        vj.input_variables.variables['ISIF'] = 3
                    else:
                        vj.input_variables.variables['ISIF'] = 3
                else:
                    vj.input_variables.variables['ISIF'] = 3
            else:
                vj.input_variables.variables['ISIF'] = 2

            # How to change IBRION
            # if info['avg_force'] < 0.1 and info['avg_stress_diag'] < 0.1 and info['avg_stress_non_diag'] < 0.1:
            #    vj.input_variables.variables['IBRION'] = 1
            # elif info['avg_force'] < 1 and info['avg_stress_diag'] < 1 and info['avg_stress_non_diag'] < 1:
            #    vj.input_variables.variables['IBRION'] = 2
            # else:
            #    vj.input_variables.variables['IBRION'] = 3

            # How to change EDIFF
            if vj.input_variables.variables['EDIFF'] > 2 * 1E-4:
                vj.input_variables.variables['EDIFF'] = round_small(vj.input_variables.variables['EDIFF'] / 2)
            else:
                vj.input_variables.variables['EDIFF'] = 1E-4

            # How to change EDIFFG
            if vj.input_variables.variables['EDIFFG'] < - 2 * self.target_force:
                vj.input_variables.variables['EDIFFG'] = round_small(vj.input_variables.variables['EDIFFG'] / 2)
            else:
                vj.input_variables.variables['EDIFFG'] = - self.target_force

            # Print new values
            print '[' + str(entry_id) + ']' + 'New Values:'
            for i in ['ISIF', 'IBRION', 'EDIFF', 'EDIFFG']:
                print '[' + str(entry_id) + ']' + i + ' : ', vj.input_variables.variables[i]
            print '-'

            for i in ['OUTCAR']:
                if not os.path.exists(workdir + os.sep + i):
                    wf = open(workdir + os.sep + i, 'w')
                    wf.write('')
                    wf.close()
                log = logging.handlers.RotatingFileHandler(workdir + os.sep + i, maxBytes=1, backupCount=1000)
                log.doRollover()

            try:
                vj.structure = read_poscar(workdir + os.sep + 'CONTCAR')
            except ValueError:
                print 'Error reading CONTCAR'

            vj.set_inputs()
            properties = vj.outcar
            status = self.status[entry_id]
            newentry = self.population.db.update(entry_id, structure=vj.structure, properties=properties, status=status)

            vj.save_json(workdir + os.sep + 'vaspjob.json')
            wf = open(workdir + os.sep + 'entry.json', 'w')
            json.dump(newentry, wf, sort_keys=True, indent=4, separators=(',', ': '))
            wf.close()
            return True
        else:
            vj.set_inputs()
            status = self.status[entry_id]
            newentry = self.population.db.update(entry_id, structure=vj.structure, status=status)

            vj.save_json(workdir + os.sep + 'vaspjob.json')
            wf = open(workdir + os.sep + 'entry.json', 'w')
            json.dump(newentry, wf, sort_keys=True, indent=4, separators=(',', ': '))
            wf.close()
            return True

    def update_history(self, entry_id):
        filename = 'pychemia_relaxation.json'
        filepath = self.basedir + os.sep + entry_id + os.sep + filename
        if not os.path.exists(filepath):
            wf = open(filepath, 'w')
            data = [self.vasp_jobs[entry_id].to_dict]
            json.dump(data, wf, sort_keys=True, indent=4, separators=(',', ': '))
            wf.close()
        else:
            rf = open(filepath, 'r')
            data = json.load(rf)
            rf.close()
            data.append(self.vasp_jobs[entry_id].to_dict)
            wf = open(filepath, 'w')
            json.dump(data, wf, sort_keys=True, indent=4, separators=(',', ': '))
            wf.close()

    @property
    def workdirs(self):
        return [self.basedir + os.sep + name for name in self.population.members]

    @property
    def active_workdirs(self):
        return [self.basedir + os.sep + name for name in self.population.actives]

    def run(self, runner):

        entries_ids = self.population.members

        def worker(workdir):
            wf = open(workdir + os.sep + 'LOCK', 'w')
            wf.write('')
            wf.close()
            runner.run()
            os.remove(workdir + os.sep + 'LOCK')

        def checker(workdir):
            if os.path.isfile(workdir + os.sep + 'LOCK'):
                return False
            return self.update(workdir)

        workdirs = [self.basedir + os.sep + i for i in self.population.actives]
        runner.run_multidirs(workdirs, worker, checker)

        if not self.is_running:
            self.run(runner)

    def set_run(self, code, runner, basedir, density_of_kpoints=10000, encut=1.1):

        self.runner = runner

        self.create_dirs(clean=True)
        self.create_inputs(density_of_kpoints=density_of_kpoints, encut=encut)


class Polarization:
    def __init__(self, structure, path, potcar_filepath, external=None, maxfield=0, stepfield=0):

        self.structure = structure
        self.path = path
        self.external = None
        self.potcar = potcar_filepath

        if external is not None and extenal in ['electric', 'magnetic']:
            self.external = external

        self.maxfield = abs(maxfield)
        self.stepfield = abs(stepfield)

    def initialize(self, kpoints, cleandir=False):
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
                pathname = self.path + os.sep + 'Field' + sign + str(field)
                if not os.path.isdir(pathname):
                    os.mkdir(pathname)
                for step in ['SCF', 'BANDS', 'berry_IGPAR1', 'berry_IGPAR2', 'berry_IGPAR3']:
                    if not os.path.isdir(pathname + os.sep + step):
                        os.mkdir(pathname + os.sep + step)
                    save_KPOINTS(kpoints, pathname + os.sep + step + os.sep + 'KPOINTS')
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
                pathname = self.path + os.sep + 'Field' + sign + str(field)
                for step in ['SCF', 'BANDS', 'berry_IGPAR1', 'berry_IGPAR2', 'berry_IGPAR3']:
                    tk = Tasks()
                    tk.minimum(ENCUT=200, POTCAR=self.potcar)
                    tk.vaspinput['EDIFF'] = 1E-9
                    tk.vaspinput['IBRION'] = -1
                    tk.vaspinput['NSW'] = 0
                    if step == 'SCF':
                        iv = InputVariables(variables=tk.vaspinput)
                        save_INCAR(iv, pathname + os.sep + step + os.sep + 'INCAR')
                        save_POSCAR(structure=self.structure, filepath=pathname + os.sep + step + os.sep + 'POSCAR')
                        options = {'nproc': 4, 'code_bin': '/home/guilleaf/local/src/vasp.5.3/vasp', 'mpi': True}
                        runner = Runner('vasp', 'local', options)
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
                        save_INCAR(iv, pathname + os.sep + step + os.sep + 'INCAR')
                        save_POSCAR(structure=self.structure, filepath=pathname + os.sep + step + os.sep + 'POSCAR')
                        options = {'nproc': 4, 'code_bin': '/home/guilleaf/local/src/vasp.5.3/vasp', 'mpi': True}
                        runner = Runner('vasp', 'local', options)
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
                        save_INCAR(iv, pathname + os.sep + step + os.sep + 'INCAR')
                        save_POSCAR(structure=self.structure, filepath=pathname + os.sep + step + os.sep + 'POSCAR')
                        options = {'nproc': 4, 'code_bin': '/home/guilleaf/local/src/vasp.5.3/vasp', 'mpi': True}
                        runner = Runner('vasp', 'local', options)
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


class Convergence:

    def __init__(self, energy_tolerance):

        self.convergence_info = None
        self.energy_tolerance = energy_tolerance

    def _best_value(self, variable):
        if not self.is_converge:
            print 'Convergence not completed'
            return None
        else:
            return self.convergence_info[-3][variable]

    @property
    def is_converge(self):
        if self.convergence_info is None or len(self.convergence_info) < 3:
            return False

        energies = [x['free_energy'] for x in self.convergence_info]
        if len(energies) > 2 and abs(max(energies[-3:])-min(energies[-3:])) >= self.energy_tolerance:
            return False
        else:
            return True

    def _convergence_plot(self, variable, xlabel, title, figname, annotate):

        if not self.is_converge:
            print 'Convergence not executed'
            return

        import matplotlib
        matplotlib.use('agg')
        import matplotlib.pyplot as plt
        x = [idata[variable] for idata in self.convergence_info]
        y = [idata['free_energy'] for idata in self.convergence_info]
        plt.clf()
        plt.plot(x, y, 'rd-')
        dy = self.energy_tolerance
        sup_dy = min(y[-3:])+dy
        low_dy = max(y[-3:])-dy
        xlims = plt.xlim()
        plt.plot(xlims, [sup_dy, sup_dy], '0.5')
        plt.plot(xlims, [low_dy, low_dy], '0.5')
        plt.fill_between(xlims, [low_dy, low_dy], [sup_dy, sup_dy], color='0.9', alpha=0.5)
        for idata in self.convergence_info:
            plt.annotate(s=str(idata[annotate]), xy=(idata[variable], idata['free_energy']), size=10)

        plt.xlim(*xlims)
        plt.title(title)
        plt.xlabel(xlabel)
        plt.ylabel('Free Energy [eV]')
        plt.savefig(figname)
        return plt.gca()

    def _convergence_save(self, filename):

        if not self.is_converge:
            print 'Convergence not executed'
            return

        wf = open(filename, 'w')
        json.dump(self.convergence_info, wf, sort_keys=True, indent=4, separators=(',', ': '))
        wf.close()

    def _convergence_load(self, filename):

        if not os.path.isfile(filename):
            raise ValueError('File not found: %s', filename)

        rf = open(filename, 'r')
        self.convergence_info = json.load(rf)
        rf.close()


class ConvergenceCutOffEnergy(Convergence):

    def __init__(self, structure, workdir='.', kpoints=None, binary='vasp', nparal=4, energy_tolerance=1E-3,
                 increment_factor=0.2, initial_encut=1.3):

        self.structure = structure
        self.workdir = workdir
        self.binary = binary
        self.nparal = nparal
        self.increment_factor = increment_factor
        self.initial_encut = initial_encut
        if kpoints is None:
            kp = KPoints()
            kp.set_optimized_grid(self.structure.lattice, density_of_kpoints=1000, force_odd=True)
            self.kpoints = kp
        else:
            self.kpoints = kpoints
        Convergence.__init__(self, energy_tolerance)

    def run(self):

        vj = VaspJob()
        vj.initialize(self.workdir, self.structure, self.kpoints, binary=self.binary)
        energies = []
        if not self.is_converge:
            x = self.initial_encut
        else:
            x = self.convergence_info[-3]['factor']
        self.convergence_info = []
        while True:
            vj.clean()
            vj.job_static()
            vj.input_variables.set_density_for_restart()
            vj.input_variables.set_encut(ENCUT=x, POTCAR=self.workdir+os.sep+'POTCAR')
            vj.set_inputs()
            encut = vj.input_variables.variables['ENCUT']
            pcm_log.debug('Testing ENCUT = %7.3f' % encut)
            vj.run(use_mpi=True, mpi_num_procs=self.nparal)
            pcm_log.debug('Starting VASP')
            while True:
                energy_str = ''
                filename = self.workdir + os.sep + 'vasp_stdout.log'
                if os.path.exists(filename):
                    vasp_stdout = read_vasp_stdout(filename=filename)
                    if len(vasp_stdout['data']) > 2:
                        scf_energies = [i[2] for i in vasp_stdout['data']]
                        energy_str = ' %7.3f' % scf_energies[1]
                        for i in range(1, len(scf_energies)):
                            if scf_energies[i] < scf_energies[i-1]:
                                energy_str += ' >'
                            else:
                                energy_str += ' <'
                        pcm_log.debug(energy_str)

                if vj.runner is not None and vj.runner.poll() is not None:
                    filename = self.workdir + os.sep + 'vasp_stdout.log'
                    if os.path.exists(filename):
                        vasp_stdout = read_vasp_stdout(filename=filename)
                        if len(vasp_stdout['data']) > 2:
                            scf_energies = [i[2] for i in vasp_stdout['data']]
                            energy_str += ' %7.3f' % scf_energies[-1]
                            pcm_log.debug(energy_str)
                    pcm_log.debug('Execution complete')
                    break
                time.sleep(5)
            vj.get_outputs()
            free_energy = vj.outcar.final_data['energy']['free_energy']
            print 'encut= %7.3f  free_energy: %9.6f' % (encut, free_energy)
            self.convergence_info.append({'free_energy': free_energy, 'encut': encut, 'factor': x})
            energies.append(free_energy)
            if len(energies) > 2 and abs(max(energies[-3:])-min(energies[-3:])) < self.energy_tolerance:
                break
            x = round(x + x*self.increment_factor, 2)

    @property
    def best_encut(self):
        return self._best_value('encut')

    def plot(self, figname='convergence_encut.pdf'):
        return self._convergence_plot(variable='encut', xlabel='ENCUT', title='ENCUT Convergence', figname=figname,
                                      annotate='encut')

    def save(self, filename='convergence_encut.json'):
        self._convergence_save(filename=filename)

    def load(self, filename='convergence_encut.json'):
        self._convergence_load(filename=filename)


class ConvergenceKPointGrid(Convergence):

    def __init__(self, structure, workdir='.', binary='vasp', nparal=4, energy_tolerance=1E-3, recover=False,
                 encut=1.3):

        self.structure = structure
        self.workdir = workdir
        self.binary = binary
        self.nparal = nparal
        self.initial_number = 12
        self.convergence_info = None
        self.encut = encut
        Convergence.__init__(self, energy_tolerance)
        if recover:
            self.recover()

    def recover(self):
        kpoints_file = self.workdir+os.sep+'KPOINTS'
        poscar_file = self.workdir+os.sep+'POSCAR'
        if os.path.isfile(kpoints_file) and os.path.isfile(poscar_file):
            structure = read_poscar(poscar_file)
            kpoints = read_kpoints(kpoints_file)
            density = kpoints.get_density_of_kpoints(structure.lattice)
            self.initial_number = int(density ** (1.0/3.0)) - 1

    def run(self):

        vj = VaspJob()
        kp = KPoints()
        vj.initialize(self.workdir, self.structure, kp, binary=self.binary)
        grid = None
        energies = []
        if not self.is_converge:
            n = self.initial_number
        else:
            n = self.convergence_info[-3]['kp_n']
        self.convergence_info = []
        while True:
            density = n**3
            kp.set_optimized_grid(self.structure.lattice, density_of_kpoints=density, force_odd=True)
            pcm_log.debug('Trial density: %d  Grid: %s' % (density, kp.grid))
            if np.sum(grid) != np.sum(kp.grid):
                grid = kp.grid
                vj.set_kpoints(kp)
                vj.clean()
                vj.job_static()
                vj.input_variables.set_density_for_restart()
                vj.input_variables.set_encut(ENCUT=self.encut, POTCAR=self.workdir+os.sep+'POTCAR')
                vj.set_inputs()
                vj.run(use_mpi=True, mpi_num_procs=self.nparal)
                while True:
                    energy_str = ''
                    filename = self.workdir + os.sep + 'vasp_stdout.log'
                    if os.path.exists(filename):
                        vasp_stdout = read_vasp_stdout(filename=filename)
                        if len(vasp_stdout['data']) > 2:
                            scf_energies = [i[2] for i in vasp_stdout['data']]
                            energy_str = ' %7.3f' % scf_energies[1]
                            for i in range(1, len(scf_energies)):
                                if scf_energies[i] < scf_energies[i-1]:
                                    energy_str += ' >'
                                else:
                                    energy_str += ' <'
                            pcm_log.debug(energy_str)

                    if vj.runner is not None and vj.runner.poll() is not None:
                        filename = self.workdir + os.sep + 'vasp_stdout.log'
                        if os.path.exists(filename):
                            vasp_stdout = read_vasp_stdout(filename=filename)
                            if len(vasp_stdout['data']) > 2:
                                scf_energies = [i[2] for i in vasp_stdout['data']]
                                energy_str += ' %7.3f' % scf_energies[-1]
                                pcm_log.debug(energy_str)
                        break
                    time.sleep(5)
                vj.get_outputs()
                energy = vj.outcar.final_data['energy']['free_energy']
                energies.append(energy)
                print 'kp_density= %10d kp_grid= %15s free_energy= %9.6f' % (density, grid, energy)
                self.convergence_info.append({'free_energy': vj.outcar.final_data['energy']['free_energy'],
                                              'kp_grid': list(grid),
                                              'kp_density': density,
                                              'kp_n': n})
                if len(energies) > 2 and abs(max(energies[-3:])-min(energies[-3:])) < self.energy_tolerance:
                    break
            n += 2

    def plot(self, figname='convergence_kpoints.pdf'):
        return self._convergence_plot(variable='kp_density', xlabel='K-points density', title='KPOINTS Convergence',
                                      figname=figname, annotate='kp_grid')

    def save(self, filename='convergence_kpoints.json'):
        self._convergence_save(filename=filename)

    def load(self, filename='convergence_kpoints.json'):
        self._convergence_load(filename=filename)

    @property
    def best_kpoints(self):

        if not self.is_converge:
            print 'Convergence not completed'
            return None
        else:
            kp = KPoints()
            kp.set_optimized_grid(self.structure.lattice, density_of_kpoints=self.convergence_info[-3]['kp_density'],
                                  force_odd=True)
            return kp


class VaspRelaxator(Relaxator):

    def __init__(self, workdir, structure, relaxator_params, target_forces=1E-3, waiting=False, binary='vasp'):

        Relaxator.__init__(self, target_forces)
        self.encut = 1.3
        self.workdir = workdir
        self.initial_structure = structure
        self.structure = self.initial_structure.copy()
        self.kpoints = None
        self.target_forces = target_forces
        self.waiting = waiting
        self.vaspjob = VaspJob()
        self.relaxed = False
        self.binary = binary
        self.nparal = None
        self.set_params(relaxator_params)
        self.vaspjob.initialize(workdir=self.workdir, structure=self.structure, kpoints=self.kpoints,
                                binary=self.binary)

    def create_dirs(self, clean=False):
        if not os.path.isdir(self.workdir):
            os.makedirs(self.workdir)
        elif clean:
            for i in os.listdir(self.workdir):
                shutil.rmtree(self.workdir + os.sep + i)

    def set_params(self, params):
        if 'kp_grid' in params:
            self.kpoints = KPoints(kmode='gamma', grid=params['kp_grid'])
        elif 'kp_density' in params:
            self.kpoints = KPoints()
            self.kpoints.set_optimized_grid(self.structure.lattice, density_of_kpoints=params['kp_density'])
        else:
            self.kpoints = KPoints()
            self.kpoints.set_optimized_grid(self.structure.lattice, density_of_kpoints=10000)
        if 'encut' in params:
            self.encut = params['encut']
        else:
            pass
        if 'nparal' in params:
            self.nparal = params['nparal']
        else:
            self.nparal = 4

    def add_status(self, status):
        pass

    def del_status(self, status):
        pass

    def update(self):
        """
        This routine determines how to proceed with the relaxation
        for one specific work directory

        :return:
        """
        vj = self.vaspjob

        if os.path.isfile(self.workdir + os.sep + 'OUTCAR'):
            vj.get_outputs()

        max_force, max_stress = self.get_max_force_stress()
        if max_force is not None and max_stress is not None:
            pcm_log.debug('Max Force: %9.3E Stress: %9.3E' % (max_force, max_stress))
            vo = VaspOutput(self.workdir + os.sep + 'OUTCAR')
            info = vo.relaxation_info()
            pcm_log.debug('Avg Force: %9.3E Stress: %9.3E %9.3E' % (info['avg_force'],
                                                                    info['avg_stress_diag'],
                                                                    info['avg_stress_non_diag']))
        else:
            print 'Failure to get forces and stress'
            return False

        if os.path.isfile(self.workdir + os.sep + 'RELAXED'):
            self.add_status('RELAXED')
        elif not os.path.isfile(self.workdir + os.sep + 'PROCAR'):
            self.add_status('NOPROCAR')
        else:
            self.del_status('NOPROCAR')
            if not os.path.isfile(self.workdir + os.sep + 'OUTCAR'):
                self.add_status('NOOUTCAR')
            else:
                self.del_status('NOOUTCAR')
                vo = VaspOutput(self.workdir + os.sep + 'OUTCAR')
                info = vo.relaxation_info()
                if len(info) != 3:
                    print ' Missing some data in OUTCAR (forces or stress)'
                    self.add_status('NOOUTCAR')

                # Conditions to consider the structure relaxed
                if info['avg_force'] < self.target_forces:
                    if info['avg_stress_diag'] < self.target_forces:
                        if info['avg_stress_non_diag'] < self.target_forces:
                            wf = open(self.workdir + os.sep + 'RELAXED', 'w')
                            for i in info:
                                wf.write("%15s %12.3f" % (i, info[i]))
                            wf.close()
                            wf = open(self.workdir + os.sep + 'COMPLETE', 'w')
                            for i in info:
                                wf.write("%15s %12.3f" % (i, info[i]))
                            wf.close()
                            self.add_status('RELAXED')

        # How to change ISIF
        if info['avg_force'] < 0.1:
            if info['avg_stress_diag'] < 0.1:
                if info['avg_stress_non_diag'] < 0.1:
                    vj.input_variables.variables['ISIF'] = 3
                else:
                    vj.input_variables.variables['ISIF'] = 3
            else:
                vj.input_variables.variables['ISIF'] = 3
        else:
            vj.input_variables.variables['ISIF'] = 2

        # How to change IBRION
        # if info['avg_force'] < 0.1 and info['avg_stress_diag'] < 0.1 and info['avg_stress_non_diag'] < 0.1:
        #    vj.input_variables.variables['IBRION'] = 1
        # elif info['avg_force'] < 1 and info['avg_stress_diag'] < 1 and info['avg_stress_non_diag'] < 1:
        #    vj.input_variables.variables['IBRION'] = 2
        # else:
        #    vj.input_variables.variables['IBRION'] = 3

        # if vj.input_variables.variables['EDIFFG'] < - 2 * self.target_forces:
        #     vj.input_variables.variables['EDIFFG'] = round_small(vj.input_variables.variables['EDIFFG'] / 2)
        # else:
        #     vj.input_variables.variables['EDIFFG'] = - self.target_forces
        #

        # How to change EDIFFG
        if max_force > self.target_forces or max_stress > self.target_forces:
            vj.input_variables.variables['EDIFFG'] = round_small(-0.01*max(max_force, max_stress))

        pcm_log.debug('Current Values: ISIF: %2d   IBRION: %2d   EDIFF: %7.1E \tEDIFFG: %7.1E' %
                      (vj.input_variables.variables['ISIF'],
                       vj.input_variables.variables['IBRION'],
                       vj.input_variables.variables['EDIFF'],
                       vj.input_variables.variables['EDIFFG']))

        # How to change EDIFF
        if vj.input_variables.variables['EDIFF'] > -0.01 * vj.input_variables.variables['EDIFFG']:
            vj.input_variables.variables['EDIFF'] = round_small(-0.01*vj.input_variables.variables['EDIFFG'])
        else:
            vj.input_variables.variables['EDIFF'] = 1E-4

        # Print new values
        pcm_log.debug('New Values: ISIF: %2d   IBRION: %2d   EDIFF: %7.1E \tEDIFFG: %7.1E' %
                      (vj.input_variables.variables['ISIF'],
                       vj.input_variables.variables['IBRION'],
                       vj.input_variables.variables['EDIFF'],
                       vj.input_variables.variables['EDIFFG']))

        for i in ['POSCAR', 'INCAR', 'OUTCAR', 'vasprun.xml']:
            if not os.path.exists(self.workdir + os.sep + i):
                wf = open(self.workdir + os.sep + i, 'w')
                wf.write('')
                wf.close()
            log = logging.handlers.RotatingFileHandler(self.workdir + os.sep + i, maxBytes=1, backupCount=1000)
            log.doRollover()

        vj.structure = self.get_final_geometry()

        vj.set_inputs()
        vj.save_json(self.workdir + os.sep + 'vaspjob.json')
        return True

    def first_run(self):

        vj = self.vaspjob
        vj.clean()
        inp = InputVariables()
        inp.set_rough_relaxation()
        vj.set_input_variables(inp)
        vj.write_potcar()
        vj.input_variables.set_encut(ENCUT=self.encut, POTCAR=self.workdir + os.sep + 'POTCAR')
        vj.input_variables.set_density_for_restart()
        vj.set_kpoints(self.kpoints)
        vj.set_inputs()

        vj.run(use_mpi=True, mpi_num_procs=self.nparal)
        if self.waiting:
            vj.runner.wait()

    def run(self):

        vj = self.vaspjob
        self.first_run()

        while True:
            if vj.runner is not None and vj.runner.poll() is not None:
                pcm_log.info('Execution completed. Return code %d' % vj.runner.returncode)

                filename = self.workdir + os.sep + 'vasp_stdout.log'
                if os.path.exists(filename):
                    read_vasp_stdout(filename=filename)

                va = VaspAnalyser(self.workdir)
                va.run()

                max_force, max_stress = self.get_max_force_stress()
                if max_force is not None and max_force < self.target_forces and max_stress < self.target_forces:
                    print 'Max Force: %9.3E Stress: %9.3E' % (max_force, max_stress)
                    break

                self.update()

                vj.run(use_mpi=True, mpi_num_procs=self.nparal)
                if self.waiting:
                    vj.runner.wait()
            else:
                filename = self.workdir + os.sep + 'vasp_stdout.log'
                if os.path.exists(filename):
                    vasp_stdout = read_vasp_stdout(filename=filename)
                    if len(vasp_stdout['iterations']) > 0:
                        pcm_log.debug('[%s] SCF: %s' % (os.path.basename(self.workdir), str(vasp_stdout['iterations'])))
                    # if len(vasp_stdout['energies']) > 2:
                    #     energy_str = ' %9.3E' % vasp_stdout['energies'][0]
                    #     for i in range(1, len(vasp_stdout['energies'])):
                    #         if vasp_stdout['energies'][i] < vasp_stdout['energies'][i-1]:
                    #             energy_str += ' >'
                    #         else:
                    #             energy_str += ' <'
                    #     pcm_log.debug(energy_str)

                time.sleep(10)

    def get_forces_stress_energy(self):

        filename = self.workdir + os.sep + 'OUTCAR'
        if os.path.isfile(filename):
            self.vaspjob.get_outputs()
            if self.vaspjob.outcar.has_forces_stress_energy():
                forces = self.vaspjob.outcar.forces[-1]
                stress = self.vaspjob.outcar.stress[-1]
                total_energy = self.vaspjob.outcar.final_data['energy']['free_energy']
            else:
                forces = None
                stress = None
                total_energy = None
        else:
            forces = None
            stress = None
            total_energy = None
        return forces, stress, total_energy

    def get_final_geometry(self):

        filename = self.workdir + os.sep + 'CONTCAR'
        structure = None
        if os.path.isfile(filename):
            try:
                structure = read_poscar(filename)
            except ValueError:
                print 'Error reading CONTCAR'
        return structure
