import os
import json
import numpy as np
from pychemia.crystal import KPoints
from ...tasks import Task
from ..abinit import AbinitJob

__author__ = 'Guillermo Avendano-Franco'


class StaticCalculation(Task):
    def __init__(self, structure, workdir='.', executable='abinit', ecut=50, kpoints=None, kp_density=1E4):

        self.ecut = ecut
        if kpoints is None:
            kp = KPoints.optimized_grid(structure.lattice, kp_density=kp_density, force_odd=True)
            self.kpoints = kp
        else:
            self.kpoints = kpoints
        self.task_params = {'ecut': self.ecut, 'kpoints': self.kpoints.to_dict}
        Task.__init__(self, structure=structure, task_params=self.task_params, workdir=workdir, executable=executable)
        self.abinitjob = AbinitJob()
        self.abinitjob.initialize(structure=structure)
        self.is_prepared=False

    def prepare(self):
        self.abinitjob.set_kpoints(kpoints=self.kpoints)
        self.abinitjob.job_static()
        self.abinitjob.set_ecut(self.ecut)
        self.abinitjob.set_psps()
        self.abinitjob.write_all()
        self.is_prepared = True

    def run(self, num_threads=None, mpi_num_procs=None):
        if not self.is_prepared:
            self.prepare()
        self.abinitjob.run(num_threads=num_threads, mpi_num_procs=mpi_num_procs)

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
        self.ecut = self.task_params['ecut']
        self.kpoints = KPoints.from_dict(self.task_params['kpoints'])

    def report(self, file_format='html'):
        from lxml.builder import ElementMaker, E
        self.plot(figname=self.report_dir + os.sep + 'static.jpg')

        element_maker = ElementMaker(namespace=None, nsmap={None: "http://www.w3.org/1999/xhtml"})
        html = element_maker.html(E.head(E.title("ABINIT Static Calculation")),
                                  E.body(E.h1("ABINIT Static Calculation"),
                                         E.h2('Structure'),
                                         E.pre(str(self.structure)),
                                         E.h2('Self Consistent Field Convergence'),
                                         E.p(E.img(src='static.jpg', width="800", height="600",
                                                   alt="Static Calculation"))
                                         ))

        return self.report_end(html, file_format)
