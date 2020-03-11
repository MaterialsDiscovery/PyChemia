import os
import json
import numpy as np
from ..abinit import AbinitJob
from ...tasks import Task
from ...relaxator import Relaxator
from pychemia.crystal import KPoints
from pychemia import pcm_log

__author__ = 'Guillermo Avendano-Franco'


class IonRelaxation(Relaxator, Task):
    def __init__(self, structure, workdir='.', tolmxf=1E-4, tolrff=1E-2, waiting=False, executable='abinit',
                 ecut=50, kp_grid=None, kp_density=1E4, relax_cell=True, max_calls=10):

        Relaxator.__init__(self, tolmxf)
        self.tolmxf = tolmxf
        self.ecut = ecut
        self.waiting = waiting
        self.abinitjob = AbinitJob(workdir=workdir, executable=executable)
        self.relaxed = False
        self.tolmxf = tolmxf
        self.tolrff = tolrff
        if kp_grid is not None:
            self.kpoints = KPoints(kmode='gamma', grid=kp_grid)
        else:
            self.kpoints = KPoints.optimized_grid(structure.lattice, kp_density=kp_density)
        self.abinitjob.initialize(structure=structure)
        self.relax_cell = relax_cell
        self.max_calls = max_calls
        task_params = {'tolmxf': self.target_forces, 'ecut': self.ecut, 'relax_cell': self.relax_cell,
                       'max_calls': self.max_calls}
        Task.__init__(self, structure=structure, task_params=task_params, workdir=workdir, executable=executable)

    def run(self, nparal=1, verbose=False):
        
        self.abinitjob.set_kpoints(kpoints=self.kpoints)
        self.abinitjob.job_ion_relax(tolmxf=self.tolmxf, tolrff=self.tolrff)
        self.abinitjob.set_ecut(self.ecut)
        self.abinitjob.set_psps()
        self.abinitjob.write_all()
        pcm_log.debug("Starting abinit execution with nparal=%d" % nparal)
        self.abinitjob.run(num_threads=nparal, mpi_num_procs=nparal, verbose=verbose)

    def get_final_geometry(self):
        pass

    def get_forces_stress_energy(self):
        pass

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
        html = element_maker.html(E.head(E.title("ABINIT Ion Relaxation")),
                                  E.body(E.h1("ABINIT Ion Relaxation"),
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
        self.ecut = self.task_params['ecut']
        self.target_forces = self.task_params['target_forces']
        self.relax_cell = self.task_params['relax_cell']
        self.max_calls = self.task_params['max_calls']
