import json
import os
import time
from pychemia import pcm_log
from pychemia.crystal import KPoints
from ..dftb import DFTBplus, read_detailed_out


class StaticCalculation:
    def __init__(self, structure, workdir, slater_path, waiting=False, kpoints=None, output_file='results.json',
                 max_scc_iterations=50):

        self.structure = structure
        self.workdir = workdir
        self.slater_path = slater_path
        self.waiting = waiting
        self.MaxSCCIterations = max_scc_iterations
        if isinstance(slater_path, str):
            self.slater_path = [slater_path]
        self.results = []
        self.output_file = output_file
        if kpoints is None:
            self.kpoints = KPoints.optimized_grid(self.structure.lattice, kp_density=10000, force_odd=True)
        else:
            self.kpoints = kpoints

    def run(self):

        dftb = DFTBplus(workdir=self.workdir)
        dftb.initialize(structure=self.structure, kpoints=self.kpoints)
        dftb.set_slater_koster(search_paths=self.slater_path)
        dftb.kpoints = self.kpoints

        if os.path.isfile('charges.bin'):
            os.remove('charges.bin')

        for mixer in ['Broyden', 'Anderson', 'DIIS', 'Simple']:

            dftb.basic_input()
            dftb.hamiltonian['MaxSCCIterations'] = self.MaxSCCIterations
            dftb.hamiltonian['Mixer'] = {'name': mixer}
            if os.path.isfile('charges.bin'):
                dftb.hamiltonian['ReadInitialCharges'] = True

            ret = None
            dftb.set_static()
            dftb.set_inputs()
            dftb.run()
            if self.waiting:
                dftb.runner.wait()
            while True:
                if dftb.runner is not None and dftb.runner.poll() is not None:
                    pcm_log.info('Execution completed. Return code %d' % dftb.runner.returncode)
                    filename = dftb.workdir + os.sep + 'detailed.out'
                    ret = read_detailed_out(filename)
                    print('Mixer= %10s  Total_energy= %9.3f  iSCC= %4d  SCC_error= %9.3E' % (mixer,
                                                                                             ret['total_energy'],
                                                                                             ret['SCC']['iSCC'],
                                                                                             ret['SCC']['SCC_error']))
                    break
                time.sleep(10)

            if ret['SCC']['iSCC'] < self.MaxSCCIterations:
                break

            if ret is not None:
                self.results.append({'Mixer'
                                     'kp_grid': self.kpoints.grid,
                                     'iSCC': ret['SCC']['iSCC'],
                                     'Total_energy': ret['total_energy'],
                                     'SCC_error': ret['SCC']['SCC_error']})

    def save_json(self):

        wf = open(self.output_file, 'w')
        json.dump(self.results, wf)
        wf.close()
