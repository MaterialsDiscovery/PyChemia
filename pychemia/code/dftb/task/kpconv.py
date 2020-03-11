import os
import time
import json
import numpy as np
from ..dftb import DFTBplus, read_detailed_out
from pychemia.crystal import KPoints
from pychemia import pcm_log, Structure


class KPointConvergence:
    def __init__(self, structure, workdir='.', slater_path='.', waiting=False, energy_tolerance=1E-3,
                 output_file='results.json'):

        self.structure = structure
        self.workdir = workdir
        self.slater_path = slater_path
        self.waiting = waiting
        self.energy_tolerance = energy_tolerance
        if isinstance(slater_path, str):
            self.slater_path = [slater_path]
        self.results = []
        self.output_file = output_file

        dftb = DFTBplus(workdir=self.workdir)
        kpoints = KPoints.optimized_grid(self.structure.lattice, kp_density=10000, force_odd=True)
        dftb.initialize(structure=self.structure, kpoints=kpoints)
        ans = dftb.set_slater_koster(search_paths=self.slater_path)
        if not ans:
            print('Slater-Koster files not complete')

    def run(self):

        n = 10
        dftb = DFTBplus(workdir=self.workdir)
        kpoints = KPoints.optimized_grid(self.structure.lattice, kp_density=10000, force_odd=True)
        dftb.initialize(structure=self.structure, kpoints=kpoints)
        ans = dftb.set_slater_koster(search_paths=self.slater_path)
        if not ans:
            print('Slater-Koster files not complete')
            return

        grid = None
        energies = []

        while True:
            density = n ** 3
            kpoints = KPoints.optimized_grid(self.structure.lattice, kp_density=density, force_odd=True)
            if np.sum(grid) != np.sum(kpoints.grid):
                pcm_log.debug('Trial density: %d  Grid: %s' % (density, kpoints.grid))
                grid = list(kpoints.grid)
                dftb.kpoints = kpoints

                dftb.basic_input()
                dftb.hamiltonian['MaxSCCIterations'] = 50
                if os.path.isfile('charges.bin'):
                    dftb.hamiltonian['ReadInitialCharges'] = True
                dftb.hamiltonian['Mixer'] = {'name': 'DIIS'}
                dftb.set_static()
                dftb.set_inputs()
                dftb.run()
                if self.waiting:
                    dftb.runner.wait()
                while True:
                    if dftb.runner is not None and dftb.runner.poll() is not None:
                        pcm_log.info('Execution completed. Return code %d' % dftb.runner.returncode)
                        filename = dftb.workdir + os.sep + 'detailed.out'
                        if os.path.exists(filename):
                            ret = read_detailed_out(filename)
                            line = 'KPoint_grid= %15s  iSCC= %4d  Total_energy= %10.4f  SCC_error= %9.3E'
                            print(line % (grid, ret['SCC']['iSCC'], ret['total_energy'], ret['SCC']['SCC_error']))
                        else:
                            print('detailed.out could not be found, exiting...')
                            return
                        n += 2
                        energies.append(ret['total_energy'])
                        break
                    time.sleep(10)

                self.results.append({'kp_grid': grid,
                                     'iSCC': ret['SCC']['iSCC'],
                                     'Total_energy': ret['total_energy'],
                                     'SCC_error': ret['SCC']['SCC_error']})
            else:
                n += 2

            if len(energies) > 2 and abs(max(energies[-3:]) - min(energies[-3:])) < self.energy_tolerance:
                break

    def save_json(self):

        wf = open(self.output_file, 'w')
        json.dump(self.results, wf, sort_keys=True, separators=(',\n', ': '))
        wf.close()


def kpoint_convergence():
    st = Structure.load_json('structure.json')
    job = KPointConvergence(st)
    job.run()
    job.save_json()
