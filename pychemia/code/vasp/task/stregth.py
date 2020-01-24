import os
import json
import pychemia
import numpy as np
from .convergence import ConvergenceKPointGrid
from .relax import IonRelaxation
from ...tasks import Task

__author__ = 'Guillermo Avendano-Franco'


class IdealStrength(Task):
    def __init__(self, structure, workdir='.', executable='vasp', ini_factor=1.0, fin_factor=1.2, delta_factor=0.01,
                 kp=None, kp_density=1E4, expansion=(1, 1, 1), encut=1.3, target_forces=1E-3, output_file=None,
                 energy_tol=1E-3):

        self.ini_factor = ini_factor
        self.fin_factor = fin_factor
        self.delta_factor = delta_factor
        self.kp_density = kp_density
        self.encut = encut
        self.target_forces = target_forces
        self.output_file = output_file
        self.expansion = expansion
        self.energy_tol = energy_tol

        if kp is None:
            kp = pychemia.crystal.KPoints.optimized_grid(structure.lattice, kp_density=kp_density, force_odd=True)
            self.kpoints = kp
        else:
            self.kpoints = kp
            self.kp_density = kp.get_density_of_kpoints(structure.lattice)

        task_params = {'ini_factor': ini_factor, 'fin_factor': fin_factor, 'delta_factor': delta_factor,
                       'kp_density': kp_density, 'expansion': expansion, 'encut': encut,
                       'target_forces': target_forces, 'kp': kp.to_dict}
        Task.__init__(self, structure=structure, task_params=task_params, executable=executable, workdir=workdir)

        self.output = []

    def run(self, nparal=4):

        for ifactor in np.arange(self.ini_factor, self.fin_factor + 0.9 * self.delta_factor, self.delta_factor):

            lattice = self.structure.lattice
            new_lengths = (ifactor - 1.0) * np.array(self.expansion) * lattice.lengths + lattice.lengths
            newlattice_params = tuple(np.concatenate(new_lengths, lattice.angles))
            newlattice = pychemia.crystal.Lattice.from_parameters_to_cell(*newlattice_params)
            newst = pychemia.Structure(cell=newlattice.cell, symbols=self.structure.symbols,
                                       reduced=self.structure.reduced)

            tmpkp = pychemia.crystal.KPoints.optimized_grid(newst.lattice, kp_density=self.kp_density, force_odd=True)

            if ifactor < 1.0:
                print('\nCompresing cell to %7.3f percent' % (ifactor * 100))
                print('-------------------\n')
            elif ifactor > 1.0:
                print('\nExpanding cell to %7.3f percent' % (ifactor * 100))
                print('-------------------\n')
            else:
                print('\nOriginal cell')
                print('-------------\n')

            new_density = tmpkp.get_density_of_kpoints(newst.lattice)
            print('KP Density: target: %d  new lattice: %d' % (self.kp_density, new_density))
            print('KP Grid:    target: %s  new lattice: %s' % (self.kpoints.grid, tmpkp.grid))

            if np.prod(tmpkp.grid) > np.prod(self.kpoints.grid):

                print('The new cell ask for a bigger grid, doing a new convergence...\n')

                self.cleaner()
                print('\nConvergence of K-Point Grid')
                print('---------------------------\n')
                ck = ConvergenceKPointGrid(newst, workdir=self.workdir + os.sep + 'KPCONV_' + str(ifactor),
                                           executable=self.executable, encut=self.encut, energy_tolerance=self.energy_tol,
                                           recover=True)
                ck.run(nparal=nparal)
                ck.save()
                ck.report()
                self.kpoints = ck.best_kpoints
            else:
                print('The current grid is still good')

            self.cleaner()
            relax = IonRelaxation(structure=newst, workdir=self.workdir + os.sep + 'RELAX_' + str(ifactor),
                                  kp_grid=self.kpoints.grid, encut=self.encut,
                                  relax_cell=False, target_forces=self.target_forces, waiting=False,
                                  executable=self.executable)
            relax.run(nparal=nparal)

            vo = pychemia.code.vasp.VaspOutput(self.workdir + '_' + str(ifactor) + '/OUTCAR')
            relst = relax.get_final_geometry()
            symm = pychemia.symm.StructureSymmetry(relst)

            self.output.append({'factor': ifactor, 'volume': newst.volume, 'density': newst.density,
                                'grid': list(self.kpoints.grid), 'output': vo.to_dict, 'spacegroup': symm.number()})
            self.save()

    def plot(self, filedir=None, file_format='pdf'):
        if filedir is None:
            filedir = self.workdir
        import matplotlib.pyplot as plt
        plt.switch_backend('agg')

        plt.figure(figsize=(8, 6))
        plt.subplots_adjust(left=0.1, bottom=0.08, right=0.95, top=0.95, wspace=None, hspace=None)
        data = np.array([[x['output']['energy'] / self.structure.natom, x['volume'] / self.output[0]['volume']]
                         for x in self.output])

        plt.plot(data[:, 0], data[:, 1])

        plt.savefig(filedir + os.sep + 'stress.' + file_format)

    def load(self, filename=None):
        if filename is None:
            filename = self.workdir + os.sep + 'task.json'
        rf = open(filename)
        data = json.load(rf)
        rf.close()
        self.task_params = data['task_params']
        self.output = data['output']
        self.ini_factor = self.task_params['ini_factor']
        self.fin_factor = self.task_params['fin_factor']
        self.delta_factor = self.task_params['delta_factor']
        self.kp_density = self.task_params['kp_density']
        self.expansion = self.task_params['expansion']
        self.encut = self.task_params['encut']
        self.target_forces = self.task_params['target_forces']
        self.kpoints = pychemia.crystal.KPoints.from_dict(self.task_params['kp'])

    def report(self, file_format='html'):
        from lxml.builder import ElementMaker, E
        self.plot(filedir=self.report_dir, file_format='jpg')

        element_maker = ElementMaker(namespace=None, nsmap={None: "http://www.w3.org/1999/xhtml"})
        html = element_maker.html(E.head(E.title("VASP Ideal Strength")),
                                  E.body(E.h1("VASP Ideal Strength"),
                                         E.h2('Structure'),
                                         E.pre(str(self.structure)),
                                         E.h2('Ideal Strength'),
                                         E.p(E.img(src='strenth.jpg', width="800", height="600", alt="Strength")),
                                         ))

        return self.report_end(html, file_format)

    def cleaner(self):

        files = [x for x in os.listdir(self.workdir) if os.path.isfile(self.workdir + os.sep + x)]
        for i in files:
            if i[:5] in ['INCAR', 'POSCA', 'OUTCA', 'vaspr']:
                os.remove(self.workdir + os.sep + i)
