import os
import json
import time
import numpy as np
from pychemia import pcm_log, HAS_MATPLOTLIB
from pychemia.crystal import KPoints
from ..vasp import VaspJob
from ..outcar import read_vasp_stdout
from ..kpoints import read_kpoints
from ..poscar import read_poscar
from ...tasks import Task
from ..xml_output import parse_vasprun


__author__ = 'Guillermo Avendano-Franco'


class Convergence:
    def __init__(self, energy_tolerance):

        self.convergence_info = []
        self.energy_tolerance = energy_tolerance

    def _best_value(self, variable):
        if not self.is_converge:
            print('Convergence not completed')
            return None
        else:
            return self.convergence_info[-3][variable]

    @property
    def is_converge(self):
        if self.convergence_info is None or len(self.convergence_info) < 3:
            return False

        energies = [x['free_energy'] for x in self.convergence_info]
        if len(energies) > 2 and abs(max(energies[-3:]) - min(energies[-3:])) >= self.energy_tolerance:
            return False
        else:
            return True

    def _convergence_plot(self, variable, xlabel, title, figname, annotate):

        import matplotlib.pyplot as plt
        if not self.is_converge:
            print('Convergence not executed')
            return

        x = [idata[variable] for idata in self.convergence_info]
        y = [idata['free_energy'] for idata in self.convergence_info]
        plt.figure(figsize=(10, 8))
        plt.clf()
        plt.plot(x, y, 'rd-')
        dy = self.energy_tolerance
        sup_dy = min(y[-3:]) + dy
        low_dy = max(y[-3:]) - dy
        xlims = plt.xlim()
        plt.plot(xlims, [sup_dy, sup_dy], '0.5')
        plt.plot(xlims, [low_dy, low_dy], '0.5')
        plt.fill_between(xlims, [low_dy, low_dy], [sup_dy, sup_dy], color='0.9', alpha=0.5)
        assert (self.convergence_info is not None)
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
            print('Convergence not executed')
            return

        wf = open(filename, 'w')
        json.dump(self.convergence_info, wf, sort_keys=True, indent=4, separators=(',', ': '))
        wf.close()

    def _convergence_load(self, filename):

        if not os.path.isfile(filename):
            raise ValueError('File not found: %s', filename)
        rf = open(filename, 'r')

        data = json.load(rf)
        self.output = data['output']
        self.task_params = data['task_params']
        self.convergence_info = self.output['convergence']
        self.energy_tolerance = self.task_params['energy_tolerance']
        rf.close()


class ConvergenceCutOffEnergy(Task, Convergence):
    def __init__(self, structure, workdir='.', kpoints=None, executable='vasp', energy_tolerance=1E-3,
                 increment_factor=0.2, initial_encut=1.3, pspdir='potpaw_PBE', psp_options=None, extra_vars=None, heterostructure=False):
        
        self.structure = structure
        self.workdir = workdir
        self.executable = executable
        self.increment_factor = increment_factor
        self.initial_encut = initial_encut
        self.pspdir = pspdir
        if psp_options is not None:
            self.psp_options = psp_options
        else:
            self.psp_options = {}
        if extra_vars is not None:
            self.extra_vars = extra_vars
        else:
            self.extra_vars = {}
        if kpoints is None:
            kp = KPoints.optimized_grid(self.structure.lattice, kp_density=1E4, force_odd=True)
            self.kpoints = kp
        else:
            self.kpoints = kpoints

        # If heterostructure is true it will keep the repeating order found
        # in the POSCAR. 
        # Added by Uthpala on Apr 20th, 2020.
        self.heterostructure = heterostructure

        Convergence.__init__(self, energy_tolerance)
        self.task_params = {'energy_tolerance': self.energy_tolerance, 'increment_factor': self.increment_factor,
                            'initial_encut': self.initial_encut}
        Task.__init__(self, structure=structure, task_params=self.task_params, workdir=workdir, executable=executable)

    def run(self, nparal=4):

        self.started = True
        vj = VaspJob(workdir=self.workdir, executable=self.executable)
        vj.initialize(structure=self.structure, kpoints=self.kpoints, pspdir=self.pspdir, heterostructure=self.heterostructure)
        vj.potcar_setup = self.psp_options
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
            vj.input_variables.set_encut(ENCUT=x, POTCAR=self.workdir + os.sep + 'POTCAR')
            vj.input_variables.variables['NBANDS'] = int(nparal * ((30 +
                                                                    self.structure.valence_electrons()) / nparal + 1))
            vj.input_variables.set_ismear(self.kpoints)
            vj.input_variables.variables['SIGMA'] = 0.2
            vj.input_variables.variables['ISPIN'] = 2
            for i in self.extra_vars:
                vj.input_variables[i] = self.extra_vars[i]
            vj.set_inputs()
            encut = vj.input_variables.variables['ENCUT']
            print('Testing ENCUT = %7.3f' % encut)
            vj.run(mpi_num_procs=nparal)
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
                            if scf_energies[i] < scf_energies[i - 1]:
                                energy_str += ' >'
                            else:
                                energy_str += ' <'
                        pcm_log.debug(energy_str)

                if vj.runner is not None and vj.runner.poll() is not None:
                    filename = self.workdir + os.sep + 'vasp_stdout.log'
                    if os.path.exists(filename):
                        # vasp_stdout = read_vasp_stdout(filename=filename)
                        vasprun = parse_vasprun('vasprun.xml')
                        if len(vasp_stdout['data']) > 2:
                            scf_energies = [i[2] for i in vasp_stdout['data']]
                            energy_str += ' %7.3f' % scf_energies[-1]
                            pcm_log.debug(energy_str)
                    pcm_log.debug('Execution complete')
                    break
                time.sleep(5)
            vj.get_outputs()
            free_energy = vj.outcar.final_data['energy']['free_energy']/self.structure.natom
            print('encut= %7.3f  free_energy: %9.6f' % (encut, free_energy))
            self.convergence_info.append({'free_energy': free_energy, 'encut': encut, 'factor': x})
            energies.append(free_energy)
            if len(energies) > 2 and abs(max(energies[-3:]) - min(energies[-3:])) < self.energy_tolerance:
                self.success = True
                break
            x = round(x + x * self.increment_factor, 2)
        self.output = {'convergence': self.convergence_info, 'best_encut': self.best_encut}
        self.finished = True

    @property
    def best_encut(self):
        return self._best_value('encut')

    def plot(self, filedir=None, file_format='pdf'):
        if filedir is None:
            filedir = self.workdir
        figname = filedir + os.sep + 'convergence.' + file_format
        return self._convergence_plot(variable='encut', xlabel='ENCUT', title='ENCUT Convergence', figname=figname,
                                      annotate='encut')

    def load(self, filename=None):
        if filename is None:
            filename = self.workdir + os.sep + 'task.json'
        self._convergence_load(filename=filename)
        self.initial_encut = self.task_params['initial_encut']
        self.increment_factor = self.task_params['increment_factor']
        self.finished = True

    def report(self, file_format='html'):

        from lxml.builder import ElementMaker, E

        if not os.path.isdir(self.report_dir):
            os.mkdir(self.report_dir)

        self.plot(filedir=self.report_dir, file_format='jpg')

        element_maker = ElementMaker(namespace=None, nsmap={None: "http://www.w3.org/1999/xhtml"})
        html = element_maker.html(E.head(E.title("VASP Energy cut-off Convergence")),
                                  E.body(E.h1("VASP Energy cut-off Convergence"),
                                         E.h2('Structure'),
                                         E.pre(str(self.structure)),
                                         E.h2('Convergence'),
                                         E.p(E.img(src='convergence.jpg', width="800", height="600", alt="Forces")),
                                         ))

        return self.report_end(html, file_format)


class ConvergenceKPointGrid(Task, Convergence):
    def __init__(self, structure, workdir='.', executable='vasp', energy_tolerance=1E-3, recover=False, encut=1.3,
                 pspdir='potpaw_PBE', extra_vars=None, psp_options=None, heterostructure=False):

        self.structure = structure
        self.workdir = workdir
        self.executable = executable
        self.initial_number = 12
        self.convergence_info = None
        self.encut = encut
        self.pspdir = pspdir
        self.psp_options = psp_options
        if extra_vars is not None:
            self.extra_vars = extra_vars
        else:
            self.extra_vars = {}

        # If heterostructure is true it will keep the repeating order found
        # in the POSCAR. 
        # Added by Uthpala on Apr 20th, 2020.
        self.heterostructure = heterostructure

        Convergence.__init__(self, energy_tolerance)
        if recover:
            self.recover()
        self.task_params = {'energy_tolerance': self.energy_tolerance, 'encut': self.encut}
        Task.__init__(self, structure=structure, task_params=self.task_params, workdir=workdir, executable=executable)

    def recover(self):
        kpoints_file = self.workdir + os.sep + 'KPOINTS'
        poscar_file = self.workdir + os.sep + 'POSCAR'
        if os.path.isfile(kpoints_file) and os.path.isfile(poscar_file):
            structure = read_poscar(poscar_file)
            kpoints = read_kpoints(kpoints_file)
            density = kpoints.get_density_of_kpoints(structure.lattice)
            self.initial_number = int(density ** (1.0 / 3.0)) - 1

    def run(self, nparal=4):

        self.started = True

        vj = VaspJob(workdir=self.workdir, executable=self.executable)
        vj.potcar_setup = self.psp_options
        kp = KPoints()
        vj.initialize(structure=self.structure, kpoints=kp, pspdir=self.pspdir, heterostructure=self.heterostructure)
        grid = None
        energies = []
        if not self.is_converge:
            n = self.initial_number
        else:
            n = self.convergence_info[-3]['kp_n']
        self.convergence_info = []
        while True:
            density = n ** 3
            kp = KPoints.optimized_grid(self.structure.lattice, kp_density=density, force_odd=True)
            pcm_log.debug('Trial density: %d  Grid: %s' % (density, kp.grid))
            if np.sum(grid) != np.sum(kp.grid):
                grid = kp.grid
                vj.set_kpoints(kp)
                vj.clean()
                vj.job_static()
                vj.input_variables.set_density_for_restart()
                vj.input_variables.set_encut(ENCUT=self.encut, POTCAR=self.workdir + os.sep + 'POTCAR')
                vj.input_variables.variables['NBANDS'] = nparal * ((30 + self.structure.valence_electrons()) //
                                                                   nparal + 1)
                vj.input_variables.set_ismear(kp)
                vj.input_variables.variables['SIGMA'] = 0.2
                vj.input_variables.variables['ISPIN'] = 2
                for i in self.extra_vars:
                    vj.input_variables[i] = self.extra_vars[i]
                vj.set_inputs()
                vj.run(mpi_num_procs=nparal)
                while True:
                    energy_str = ''
                    filename = self.workdir + os.sep + 'vasp_stdout.log'
                    if os.path.exists(filename):
                        vasp_stdout = read_vasp_stdout(filename=filename)
                        if len(vasp_stdout['data']) > 2:
                            scf_energies = [i[2] for i in vasp_stdout['data']]
                            energy_str = ' %7.3f' % scf_energies[1]
                            for i in range(1, len(scf_energies)):
                                if scf_energies[i] < scf_energies[i - 1]:
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
                energy = vj.outcar.final_data['energy']['free_energy']/self.structure.natom
                energies.append(energy)
                print('kp_density= %10d kp_grid= %15s free_energy= %9.6f' % (density, grid, energy))
                self.convergence_info.append({'free_energy': vj.outcar.final_data['energy']['free_energy']/self.structure.natom,
                                              'kp_grid': list(grid),
                                              'kp_density': density,
                                              'kp_n': n})
                if len(energies) > 2 and abs(max(energies[-3:]) - min(energies[-3:])) < self.energy_tolerance:
                    self.success = True
                    break
            n += 2
        self.output = {'convergence': self.convergence_info, 'best_kp_grid': list(self.best_kpoints.grid)}
        self.finished = True

    def plot(self, filedir=None, file_format='pdf'):
        if filedir is None:
            filedir = self.workdir
        figname = filedir + os.sep + 'convergence.' + file_format

        return self._convergence_plot(variable='kp_density', xlabel='K-points density', title='KPOINTS Convergence',
                                      figname=figname, annotate='kp_grid')

    def load(self, filename=None):
        if filename is None:
            filename = self.workdir + os.sep + 'task.json'
        self._convergence_load(filename=filename)
        self.encut = self.task_params['encut']
        self.finished = True

    @property
    def best_kpoints(self):

        if not self.is_converge:
            print('Convergence not completed')
            return None
        else:
            kp = KPoints.optimized_grid(self.structure.lattice, kp_density=self.convergence_info[-3]['kp_density'],
                                        force_odd=True)
            return kp

    def report(self, file_format='html'):
        from lxml.builder import ElementMaker, E

        if not os.path.isdir(self.report_dir):
            os.mkdir(self.report_dir)

        self.plot(filedir=self.report_dir, file_format='jpg')

        element_maker = ElementMaker(namespace=None, nsmap={None: "http://www.w3.org/1999/xhtml"})
        html = element_maker.html(E.head(E.title("VASP K-point grid Convergence")),
                                  E.body(E.h1("VASP K-point grid Convergence"),
                                         E.h2('Structure'),
                                         E.pre(str(self.structure)),
                                         E.h2('Convergence'),
                                         E.p(E.img(src='convergence.jpg', width="800", height="600", alt="Forces")),
                                         ))

        return self.report_end(html, file_format)
