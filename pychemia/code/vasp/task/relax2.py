import json
import logging.handlers
import os
import shutil
import time

import numpy as np

from pychemia import pcm_log
from pychemia.crystal import KPoints
from pychemia.utils.serializer import generic_serializer
from ..outcar import VaspOutput, read_vasp_stdout
from ..poscar import read_poscar
from ..vasp import VaspJob, VaspAnalyser, VaspOutput
from ..input import VaspInput
from ...relaxator import Relaxator
from ...tasks import Task


__author__ = 'Guillermo Avendano-Franco'


class IonRelaxation2(Relaxator, Task):
    def __init__(self, structure, workdir='.', target_forces=1E-3, waiting=False, executable='vasp',
                 encut=1.3, kp_grid=None, kp_density=1E4, relax_cell=True):

        Relaxator.__init__(self, target_forces)
        self.target_forces = target_forces
        self.waiting = waiting
        self.vaspjob = VaspJob(workdir=workdir, executable=executable)
        self.relaxed = False
        if kp_grid is not None:
            self.kpoints = KPoints(kmode='gamma', grid=kp_grid)
        else:
            self.kpoints = KPoints.optimized_grid(structure.lattice, kp_density=kp_density)
        self.vaspjob.initialize(structure=structure, kpoints=self.kpoints)
        self.encut = encut
        self.relax_cell = relax_cell
        task_params = {'target_forces': self.target_forces, 'encut': self.encut, 'relax_cell': self.relax_cell}
        Task.__init__(self, structure=structure, task_params=task_params, workdir=workdir, executable=executable)
        self.stage = 1

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

        if self.stage == 1:

            vj.input_variables.variables['PREC'] = 'LOW'
            vj.input_variables.variables['EDIFF'] = 1e-2
            vj.input_variables.variables['EDIFFG'] = 1e-1
            vj.input_variables.set_encut(1.0, POTCAR=self.workdir + os.sep + 'POTCAR')
            vj.input_variables.variables['NSW'] = 65
            vj.input_variables.variables['ISIF'] = 4
            vj.input_variables.variables['IBRION'] = 2
            vj.input_variables.variables['POTIM'] = 0.02
            vj.input_variables.variables['ISMEAR'] = 1
            vj.input_variables.variables['SIGMA'] = 0.10

        elif self.stage == 2:
            vj.input_variables.variables['PREC'] = 'NORMAL'
            vj.input_variables.variables['EDIFF'] = 1e-3
            vj.input_variables.variables['EDIFFG'] = 1e-2
            vj.input_variables.set_encut(1.1, POTCAR=self.workdir + os.sep + 'POTCAR')
            vj.input_variables.variables['NSW'] = 55
            vj.input_variables.variables['ISIF'] = 4
            vj.input_variables.variables['IBRION'] = 1
            vj.input_variables.variables['POTIM'] = 0.30
            vj.input_variables.variables['ISMEAR'] = 1
            vj.input_variables.variables['SIGMA'] = 0.08

        elif self.stage == 3:
            vj.input_variables.variables['PREC'] = 'NORMAL'
            vj.input_variables.variables['EDIFF'] = 1e-3
            vj.input_variables.variables['EDIFFG'] = 1e-2
            vj.input_variables.set_encut(1.2, POTCAR=self.workdir + os.sep + 'POTCAR')  # Originally     ENCUT=520.0
            vj.input_variables.variables['NSW'] = 65
            vj.input_variables.variables['ISIF'] = 3
            vj.input_variables.variables['IBRION'] = 2
            vj.input_variables.variables['POTIM'] = 0.02
            vj.input_variables.variables['ISMEAR'] = 1
            vj.input_variables.variables['SIGMA'] = 0.07

        elif self.stage == 4:
            vj.input_variables.variables['PREC'] = 'NORMAL'
            vj.input_variables.variables['EDIFF'] = 1e-4
            vj.input_variables.variables['EDIFFG'] = 1e-3
            vj.input_variables.set_encut(1.3, POTCAR=self.workdir + os.sep + 'POTCAR')  # Originally     ENCUT=600.0
            vj.input_variables.variables['NSW'] = 55
            vj.input_variables.variables['ISIF'] = 3
            vj.input_variables.variables['IBRION'] = 1
            vj.input_variables.variables['POTIM'] = 0.30
            vj.input_variables.variables['ISMEAR'] = 1
            vj.input_variables.variables['SIGMA'] = 0.06

        elif self.stage == 5:
            vj.input_variables.variables['PREC'] = 'NORMAL'
            vj.input_variables.variables['EDIFF'] = 1e-4
            vj.input_variables.variables['EDIFFG'] = 1e-3
            vj.input_variables.set_encut(1.4, POTCAR=self.workdir + os.sep + 'POTCAR')  # Originally     ENCUT=600.0
            vj.input_variables.variables['NSW'] = 0
            vj.input_variables.variables['ISIF'] = 2
            vj.input_variables.variables['IBRION'] = 2
            vj.input_variables.variables['POTIM'] = 0.02
            vj.input_variables.variables['ISMEAR'] = 1
            vj.input_variables.variables['SIGMA'] = 0.05

        elif self.stage == 6:
            vj.input_variables.variables['PREC'] = 'NORMAL'
            vj.input_variables.variables['EDIFF'] = 0.1 * self.target_forces
            vj.input_variables.variables['EDIFFG'] = 0.5 * self.target_forces
            vj.input_variables.set_encut(1.5, POTCAR=self.workdir + os.sep + 'POTCAR')  # Originally     ENCUT=600.0
            vj.input_variables.variables['NSW'] = 0
            vj.input_variables.variables['ISIF'] = 2
            vj.input_variables.variables['IBRION'] = 2
            vj.input_variables.variables['POTIM'] = 0.02
            vj.input_variables.variables['ISMEAR'] = 1
            vj.input_variables.variables['SIGMA'] = 0.05

        pcm_log.debug('STAGE: %d Current Values: ISIF: %2d   IBRION: %2d   EDIFF: %7.1E \tEDIFFG: %7.1E' %
                      (self.stage,
                       vj.input_variables.variables['ISIF'],
                       vj.input_variables.variables['IBRION'],
                       vj.input_variables.variables['EDIFF'],
                       vj.input_variables.variables['EDIFFG']))

        vj.input_variables.variables['LWAVE'] = False

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

        self.stage += 1

        return True

    def first_run(self, nparal=4):
        print('FIRST RUN')
        print('=========')
        vj = self.vaspjob
        vj.clean()
        inp = VaspInput()
        inp.set_rough_relaxation()
        vj.set_input_variables(inp)
        vj.input_variables.variables['PREC'] = 'LOW'
        vj.input_variables.variables['EDIFF'] = 1e-2
        vj.input_variables.variables['EDIFFG'] = 1e-1
        vj.input_variables.variables['NSW'] = 65
        vj.input_variables.variables['ISIF'] = 4
        vj.input_variables.variables['IBRION'] = 2
        vj.input_variables.variables['POTIM'] = 0.02
        vj.input_variables.variables['ISMEAR'] = 1
        vj.input_variables.variables['SIGMA'] = 0.10
        vj.input_variables.variables['LWAVE'] = False
        vj.write_potcar()
        vj.input_variables.set_encut(ENCUT=self.encut, POTCAR=self.workdir + os.sep + 'POTCAR')
        vj.input_variables.set_density_for_restart()
        vj.set_kpoints(self.kpoints)
        vj.set_inputs()
        print('Launching VASP using %d processes' % nparal)
        vj.run(num_threads=nparal, mpi_num_procs=nparal, nodefile=None, wait=True, verbose=False)
        if self.waiting:
            vj.runner.wait()

    def cleaner(self):

        files = [x for x in os.listdir(self.workdir) if os.path.isfile(self.workdir + os.sep + x)]
        for i in files:
            if i[:5] in ['INCAR', 'POSCA', 'OUTCA', 'vaspr']:
                os.remove(self.workdir + os.sep + i)

    def run(self, nparal=1):

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

                print('READING FORCES AND STRESS')

                max_force, max_stress = self.get_max_force_stress()
                if max_force is not None and max_stress is not None:
                    pcm_log.debug('Max Force: %9.3E Stress: %9.3E' % (max_force, max_stress))
                    vo = VaspOutput(self.workdir + os.sep + 'OUTCAR')
                    info = vo.relaxation_info()
                    pcm_log.debug('Avg Force: %9.3E Stress: %9.3E %9.3E' % (info['avg_force'],
                                                                            info['avg_stress_diag'],
                                                                            info['avg_stress_non_diag']))
                    if self.stage == 7 and info['avg_force'] < self.target_forces and \
                            info['avg_stress_diag'] < self.target_forces and \
                            info['avg_stress_non_diag'] < self.target_forces:
                        break

                else:
                    print('Failure to get forces and stress')

                print('UPDATING INPUTS FOR NEXT RUN')

                self.update()

                vj.run(num_threads=nparal, mpi_num_procs=nparal, nodefile=None, wait=True, verbose=False)
                if self.waiting:
                    vj.runner.wait()

            else:
                filename = self.workdir + os.sep + 'vasp_stdout.log'
                if os.path.exists(filename):
                    vasp_stdout = read_vasp_stdout(filename=filename)
                    if len(vasp_stdout['iterations']) > 0:
                        pcm_log.debug('[%s] SCF: %s' % (os.path.basename(self.workdir), str(vasp_stdout['iterations'])))

                # if os.path.isfile(self.workdir + os.sep + 'OUTCAR'):
                #     vj.get_outputs()

                time.sleep(30)

        # outcars = sorted([x for x in os.listdir(self.workdir) if x.startswith('OUTCAR')])[::-1]
        # vo = VaspOutput(self.workdir + os.sep + outcars[0])
        # forces = vo.forces
        # stress = vo.stress
        # if len(outcars) > 1:
        #     for i in outcars[1:]:
        #         vo = VaspOutput(self.workdir + os.sep + i)
        #         forces = np.concatenate((forces, vo.forces))
        #         stress = np.concatenate((stress, vo.stress))

#        vj.get_outputs()
        # self.output = {'forces': generic_serializer(forces), 'stress': generic_serializer(stress),
        #                'energy': vj.outcar.energy, 'energies': generic_serializer(vj.outcar.energies)}
#        if vj.outcar.is_finished:
#            self.finished = True

    def get_forces_stress_energy(self):

        filename = self.workdir + os.sep + 'OUTCAR'
        if os.path.isfile(filename):
            vo = VaspOutput(filename)
            if vo.has_forces_stress_energy():
                forces = vo.forces[-1]
                stress = vo.stress[-1]
                total_energy = vo.final_data['energy']['free_energy']
            else:
                print('ERROR: VaspOuput says no forces')
                forces = None
                stress = None
                total_energy = None
        else:
            print('ERROR: No OUTCAR')
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

