import os
import time
import json
import numpy as np
from pychemia.crystal import KPoints
from pychemia import pcm_log
from pychemia.utils.serializer import generic_serializer
from ..vasp import VaspJob
from ..outcar import read_vasp_stdout
from ...tasks import Task

__author__ = 'Guillermo Avendano-Franco'


class StaticCalculation(Task):
    def __init__(self, structure, workdir='.', executable='vasp', encut=1.3, kpoints=None, kp_density=1E4,
                 extra_incar=None):

        self.encut = encut
        if kpoints is None:
            kp = KPoints.optimized_grid(structure.lattice, kp_density=kp_density, force_odd=True)
            self.kpoints = kp
        else:
            self.kpoints = kpoints
        self.task_params = {'encut': self.encut, 'kpoints': self.kpoints.to_dict, 'extra_incar': extra_incar}
        Task.__init__(self, structure=structure, task_params=self.task_params, workdir=workdir, executable=executable)

    def run(self, nparal=4):

        vj = VaspJob(workdir=self.workdir, executable=self.executable)
        vj.initialize(self.structure, self.kpoints)
        vj.clean()
        vj.job_static()
        vj.input_variables.set_density_for_restart()
        vj.input_variables.set_encut(ENCUT=self.encut, POTCAR=self.workdir + os.sep + 'POTCAR')
        vj.input_variables.variables['NBANDS'] = int(nparal * ((int(self.structure.valence_electrons()) +
                                                                self.structure.natom) / nparal + 1))
        vj.input_variables.set_ismear(self.kpoints)
        vj.input_variables.variables['SIGMA'] = 0.2
        vj.input_variables.variables['ISPIN'] = 2
        if self.task_params['extra_incar'] is not None:
            for i in self.task_params['extra_incar']:
                vj.input_variables.variables[i] = self.task_params['extra_incar'][i]
        vj.set_inputs()
        self.encut = vj.input_variables.variables['ENCUT']
        pcm_log.debug('Executing VASP at %s with %d cores' % (self.executable, nparal))
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
                pcm_log.debug('VASP Execution finished')
                break
            time.sleep(5)
        vj.get_outputs()

        self.output = {'forces': generic_serializer(vj.outcar.forces), 'stress': generic_serializer(vj.outcar.stress),
                       'energy': vj.outcar.energy, 'energies': generic_serializer(vj.outcar.energies),
                       'INCAR': vj.input_variables.variables}
        if vj.outcar.is_finished:
            self.finished = True

    def plot(self, figname='static_calculation.pdf'):
        if not self.finished:
            print('The task is not finished')
            return
        import matplotlib.pyplot as plt
        plt.switch_backend('agg')
        plt.figure(figsize=(8, 6))
        plt.subplots_adjust(left=0.09, bottom=0.08, right=0.95, top=0.95, wspace=None, hspace=None)
        data = np.array(self.output['energies'])
        plt.plot(data[:, 1], data[:, 2], 'b.-')
        plt.xlabel('SCF cycle')
        plt.ylabel('Energy [eV]')

        a = plt.axes([.6, .6, .3, .3], axisbg='0.9')
        a.semilogy(data[:, 1], data[:, 2] - np.min(data[:, 2]))
        a.set_title('min energy %7.3f eV' % np.min(data[:, 2]))

        if figname is not None:
            plt.savefig(figname)
        return plt.gcf()

    def load(self, filename=None):
        if filename is None:
            filename = self.workdir + os.sep + 'task.json'
        rf = open(filename)
        data = json.load(rf)
        rf.close()
        self.task_params = data['task_params']
        self.output = data['output']
        self.encut = self.task_params['encut']
        self.kpoints = KPoints.from_dict(self.task_params['kpoints'])

    def report(self, file_format='html'):
        from lxml.builder import ElementMaker, E
        self.plot(figname=self.report_dir + os.sep + 'static.jpg')

        element_maker = ElementMaker(namespace=None, nsmap={None: "http://www.w3.org/1999/xhtml"})
        html = element_maker.html(E.head(E.title("VASP Static Calculation")),
                                  E.body(E.h1("VASP Static Calculation"),
                                         E.h2('Structure'),
                                         E.pre(str(self.structure)),
                                         E.h2('Self Consistent Field Convergence'),
                                         E.p(E.img(src='static.jpg', width="800", height="600",
                                                   alt="Static Calculation"))
                                         ))

        return self.report_end(html, file_format)
