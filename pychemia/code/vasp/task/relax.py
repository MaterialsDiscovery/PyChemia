import json
import logging.handlers
import os
import shutil
import time

import numpy as np

from pychemia import pcm_log
from pychemia.crystal import KPoints
from pychemia.utils.serializer import generic_serializer
from pychemia.utils.mathematics import round_small
from ..input import VaspInput
from ..outcar import VaspOutput, read_vasp_stdout
from ..poscar import read_poscar
from ..vasp import VaspJob, VaspAnalyser
from ...relaxator import Relaxator
from ...tasks import Task

__author__ = 'Guillermo Avendano-Franco'


class IonRelaxation(Relaxator, Task):
    def __init__(self, structure, workdir='.', target_forces=1E-3, executable='vasp',
                 encut=1.3, kp_grid=None, kp_density=1E4, relax_cell=True,
                 max_calls=10,pspdir='potpaw_PBE', psp_options=None, extra_vars=None, heterostructure=False):

        Relaxator.__init__(self, target_forces)
        self.target_forces = target_forces

        # If heterostructure is true it will keep the repeating order found
        # in the POSCAR. 
        # Added by Uthpala on Apr 20th, 2020.
        self.heterostructure = heterostructure 

        self.vaspjob = VaspJob(executable=executable, workdir=workdir)
        self.relaxed = False
        if kp_grid is not None:
            self.kpoints = KPoints(kmode='gamma', grid=kp_grid)
        else:
            self.kpoints = KPoints.optimized_grid(structure.lattice, kp_density=kp_density)
        self.vaspjob.initialize(structure=structure, kpoints=self.kpoints, pspdir=pspdir, heterostructure=self.heterostructure)
        self.vaspjob.potcar_setup = psp_options
        self.encut = encut
        self.relax_cell = relax_cell
        self.max_calls = max_calls
        self.pspdir = pspdir
        if extra_vars is not None:
            self.extra_vars = extra_vars
        else:
            self.extra_vars = {}

        

        task_params = {'target_forces': self.target_forces, 'encut': self.encut, 'relax_cell': self.relax_cell,
                       'max_calls': self.max_calls}
        Task.__init__(self, structure=structure, task_params=task_params, workdir=workdir, executable=executable)

    def create_dirs(self, clean=False):
        if not os.path.isdir(self.workdir):
            os.makedirs(self.workdir)
        elif clean:
            for i in os.listdir(self.workdir):
                shutil.rmtree(self.workdir + os.sep + i)

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
            print('Failure to get forces and stress')
            return False

        vo = VaspOutput(self.workdir + os.sep + 'OUTCAR')
        info = vo.relaxation_info()
        if len(info) != 3:
            print(' Missing some data in OUTCAR (forces or stress)')

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

        # How to change ISIF
        if info['avg_force'] < 0.1 and self.relax_cell:
            if info['avg_stress_diag'] < 0.1:
                if info['avg_stress_non_diag'] < 0.1:
                    vj.input_variables['ISIF'] = 3
                else:
                    vj.input_variables['ISIF'] = 3
            else:
                vj.input_variables['ISIF'] = 3
        else:
            vj.input_variables['ISIF'] = 2

        # How to change IBRION
        # if info['avg_force'] < 0.1 and info['avg_stress_diag'] < 0.1 and info['avg_stress_non_diag'] < 0.1:
        #    vj.input_variables['IBRION'] = 1
        # elif info['avg_force'] < 1 and info['avg_stress_diag'] < 1 and info['avg_stress_non_diag'] < 1:
        #    vj.input_variables['IBRION'] = 2
        # else:
        #    vj.input_variables['IBRION'] = 3

        # if vj.input_variables['EDIFFG'] < - 2 * self.target_forces:
        #     vj.input_variables['EDIFFG'] = round_small(vj.input_variables['EDIFFG'] / 2)
        # else:
        #     vj.input_variables['EDIFFG'] = - self.target_forces
        #

        # How to change EDIFFG
        if max_force > self.target_forces or max_stress > self.target_forces:
            if self.relax_cell:
                vj.input_variables['EDIFFG'] = min(round_small(-0.01 * max(max_force, max_stress)),
                                                       -self.target_forces)
            else:
                vj.input_variables['EDIFFG'] = min(round_small(-0.01 * max_force), -self.target_forces)

        pcm_log.debug('Current Values: ISIF: %2d   IBRION: %2d   EDIFF: %7.1E \tEDIFFG: %7.1E' %
                      (vj.input_variables['ISIF'],
                       vj.input_variables['IBRION'],
                       vj.input_variables['EDIFF'],
                       vj.input_variables['EDIFFG']))

        # How to change EDIFF
        if vj.input_variables['EDIFF'] > -0.01 * vj.input_variables['EDIFFG']:
            vj.input_variables['EDIFF'] = round_small(-0.01 * vj.input_variables['EDIFFG'])
        else:
            vj.input_variables['EDIFF'] = 1E-4

        vj.input_variables['LWAVE'] = False

        # Print new values
        pcm_log.debug('New Values: ISIF: %2d   IBRION: %2d   EDIFF: %7.1E \tEDIFFG: %7.1E' %
                      (vj.input_variables['ISIF'],
                       vj.input_variables['IBRION'],
                       vj.input_variables['EDIFF'],
                       vj.input_variables['EDIFFG']))

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

    def first_run(self, nparal=4, waiting=True):

        vj = self.vaspjob
        vj.clean()
        inp = VaspInput()
        inp.set_rough_relaxation()
        vj.set_input_variables(inp)
        vj.input_variables['LWAVE'] = False
        vj.write_potcar()
        vj.input_variables.set_encut(ENCUT=self.encut, POTCAR=self.workdir + os.sep + 'POTCAR')
        vj.input_variables.set_density_for_restart()

        for i in self.extra_vars:
            vj.input_variables[i]=self.extra_vars[i]

        vj.set_kpoints(self.kpoints)
        vj.set_inputs()
        print('Launching VASP using %d processes' % nparal)
        vj.run(mpi_num_procs=nparal, wait=waiting)

    def cleaner(self):

        files = [x for x in os.listdir(self.workdir) if os.path.isfile(self.workdir + os.sep + x)]
        for i in files:
            if i[:5] in ['INCAR', 'POSCA', 'OUTCA', 'vaspr']:
                os.remove(self.workdir + os.sep + i)

    def run(self, nparal=1, waiting=True):

        self.started = True
        self.cleaner()
        vj = self.vaspjob
        ncalls = 1
        self.first_run(nparal)

        while True:
            if vj.runner is not None and vj.runner.poll() is not None:
                pcm_log.info('Execution completed. Return code %d' % vj.runner.returncode)

                filename = self.workdir + os.sep + 'vasp_stdout.log'
                if os.path.exists(filename):
                    read_vasp_stdout(filename=filename)

                ncalls += 1
                va = VaspAnalyser(self.workdir)
                va.run()

                max_force, max_stress = self.get_max_force_stress()
                print('Max Force: %9.3E Stress: %9.3E (target forces= %E)' %
                      (max_force, max_stress, self.target_forces))

                if max_force is not None and max_force < self.target_forces:

                    # Conditions to finish the run
                    if max_stress < self.target_forces:
                        self.success = True
                        break
                    elif not self.relax_cell:
                        self.success = True
                        break
                    elif ncalls >= self.max_calls:
                        self.success = False
                        break

                self.update()

                vj.run(mpi_num_procs=nparal)
                if waiting:
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

                time.sleep(30)

        outcars = sorted([x for x in os.listdir(self.workdir) if x.startswith('OUTCAR')])[::-1]
        vo = VaspOutput(self.workdir + os.sep + outcars[0])
        forces = vo.forces
        stress = vo.stress
        if len(outcars) > 1:
            for i in outcars[1:]:
                vo = VaspOutput(self.workdir + os.sep + i)
                forces = np.concatenate((forces, vo.forces))
                stress = np.concatenate((stress, vo.stress))

        vj.get_outputs()
        self.output = {'forces': generic_serializer(forces), 'stress': generic_serializer(stress),
                       'energy': vj.outcar.energy, 'energies': generic_serializer(vj.outcar.energies)}
        if vj.outcar.is_finished:
            self.finished = True

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
                print('Error reading CONTCAR')
        return structure

    def plot(self, filedir=None, file_format='pdf'):
        if filedir is None:
            filedir = self.workdir
        import matplotlib.pyplot as plt
        plt.switch_backend('agg')

        plt.figure(figsize=(8, 6))
        plt.subplots_adjust(left=0.1, bottom=0.08, right=0.95, top=0.95, wspace=None, hspace=None)
        forces = np.array(self.output['forces'])
        maxforce = [np.max(np.apply_along_axis(np.linalg.norm, 1, x)) for x in forces]
        avgforce = [np.mean(np.apply_along_axis(np.linalg.norm, 1, x)) for x in forces]

        if np.max(maxforce) > 0.0 and np.max(avgforce) > 0.0:
            plt.semilogy(maxforce, 'b.-', label='Max force')
            plt.semilogy(avgforce, 'r.-', label='Mean force')
        else:
            plt.plot(maxforce, 'b.-', label='Max force')
            plt.plot(avgforce, 'r.-', label='Mean force')
        plt.xlabel('Ion movement iteration')
        plt.ylabel('Max Force')
        plt.savefig(filedir + os.sep + 'forces.' + file_format)
        plt.clf()

        plt.figure(figsize=(8, 6))
        plt.subplots_adjust(left=0.1, bottom=0.08, right=0.95, top=0.95, wspace=None, hspace=None)
        stress = np.array(self.output['stress'])
        diag_stress = [np.trace(np.abs(x)) for x in stress]
        offdiag_stress = [np.sum(np.abs(np.triu(x, 1).flatten())) for x in stress]
        plt.semilogy(diag_stress, 'b.-', label='diagonal')
        plt.semilogy(offdiag_stress, 'r.-', label='off-diagonal')
        plt.legend()
        plt.xlabel('Ion movement iteration')
        plt.ylabel(r'$\sum |stress|$ (diag, off-diag)')
        plt.savefig(filedir + os.sep + 'stress.' + file_format)

    def report(self, file_format='html'):
        from lxml.builder import ElementMaker, E

        self.plot(filedir=self.report_dir, file_format='jpg')

        element_maker = ElementMaker(namespace=None, nsmap={None: "http://www.w3.org/1999/xhtml"})
        html = element_maker.html(E.head(E.title("VASP Ion Relaxation")),
                                  E.body(E.h1("VASP Ion Relaxation"),
                                         E.h2('Initial Structure'),
                                         E.pre(str(self.structure)),
                                         E.h2('Forces Minimization'),
                                         E.p(E.img(src='forces.jpg', width="800", height="600", alt="Forces")),
                                         E.h2('Stress Minimization'),
                                         E.p(E.img(src='stress.jpg', width="800", height="600", alt="Stress"))
                                         ))

        return self.report_end(html, file_format)

    def load(self, filename=None):
        if filename is None:
            filename = self.workdir + os.sep + 'task.json'
        rf = open(filename)
        data = json.load(rf)
        rf.close()
        self.task_params = data['task_params']
        self.output = data['output']
        self.encut = self.task_params['encut']
        self.target_forces = self.task_params['target_forces']
        self.relax_cell = self.task_params['relax_cell']
        self.max_calls = self.task_params['max_calls']
