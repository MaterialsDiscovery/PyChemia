__author__ = 'Guillermo Avendano-Franco'

import os
import time
import numpy as np
from pychemia import pcm_log
import pychemia

try:
    from pychemia.symm import symmetrize
except ImportError:
    def symmetrize(structure):
        return structure
from pychemia.code import Relaxator


class VaspRelaxator(Relaxator):

    def get_final_geometry(self):
        struct = pychemia.code.vasp.read_poscar(self.workdir+os.sep+'CONTCAR')
        return struct

    def __init__(self, workdir, structure, relaxator_params, target_forces=1E-3, waiting=False):

        self.workdir = workdir
        self.initial_structure = symmetrize(structure)
        self.structure = self.initial_structure.copy()
        self.set_params(relaxator_params)
        self.waiting = waiting
        Relaxator.__init__(self, target_forces)

    def set_params(self, params):
        pass

    def create_inputs(self, density_of_kpoints=10000, ENCUT=1.0):
        # kpoints = pychemia.dft.KPoints(kmode='gamma', grid=[4, 4, 4])
        kpoints = pychemia.dft.KPoints()
        kpoints.set_optimized_grid(self.structure.lattice, density_of_kpoints=density_of_kpoints)
        print kpoints
        vj = pychemia.code.vasp.VaspJob()
        vj.initialize(workdir=self.workdir, structure=self.structure, kpoints=kpoints)
        inp = pychemia.code.vasp.InputVariables()
        inp.set_rough_relaxation()
        vj.set_input_variables(inp)
        vj.write_potcar()
        vj.input_variables.set_encut(ENCUT=ENCUT, POTCAR=self.workdir + os.sep + 'POTCAR')
        vj.set_inputs()
        return vj

    def run(self):
        vj = self.create_inputs(density_of_kpoints=10000, ENCUT=1.0)
        vj.run()
        if self.waiting:
            vj.runner.wait()
        while True:
            if vj.runner is not None and vj.runner.poll() is not None:
                pcm_log.info('Execution completed. Return code %d' % vj.runner.returncode)
                filename = self.workdir + os.sep + 'OUTCAR'
                if not os.path.exists(filename):
                    pcm_log.error('Could not find ' + filename)
                    break
                forces, stress, total_energy = self.get_forces_stress_energy()

                # Forces are good if we are at least one order of magnitude higher than target
                if forces is not None and np.max(np.abs(forces.flatten())) < 10 * self.target_forces:
                    good_forces = True
                else:
                    if forces is not None:
                        pcm_log.info('Forces: ' + str(np.max(np.abs(forces.flatten()), 10 * self.target_forces)))
                    good_forces = False

                # Stress is good if we are at least one order of magnitude higher than target
                if stress is not None and np.max(np.abs(stress.flatten())) < 10 * self.target_forces:
                    good_stress = True
                else:
                    if stress is not None:
                        pcm_log.info('Stress: ' + str(np.max(np.abs(stress.flatten()), 10 * self.target_forces)))
                    good_stress = False

                score = self.quality(dftb, score)
                pcm_log.debug('Present score : ' + str(score))
                if score < 20:

                    if good_forces and good_stress:
                        pcm_log.debug('Convergence: Internals + Cell')
                        dftb.driver['MovedAtoms'] = '1:-1'
                        dftb.driver['LatticeOpt'] = True
                    elif not good_forces and good_stress:
                        pcm_log.debug('Convergence: Internals')
                        dftb.driver['LatticeOpt'] = False
                        dftb.driver['MovedAtoms'] = '1:-1'
                    elif good_forces and not good_stress:
                        pcm_log.debug('Convergence: Internals + Cell')
                        dftb.driver['LatticeOpt'] = True
                        dftb.driver['MovedAtoms'] = '1:-1'

                    dftb.structure = read_geometry_gen(dftb.workdir + os.sep + 'geo_end.gen')
                    dftb.structure = symmetrize(dftb.structure)
                    self.structure = dftb.structure

                    dftb.get_geometry()
                    dftb.roll_outputs(irun)
                    dftb.set_inputs()
                    irun += 1
                    dftb.run()
                    if self.waiting:
                        dftb.runner.wait()
                else:
                    pcm_log.debug('Final static calculation')
                    dftb.structure = self.get_final_geometry()
                    dftb.structure = symmetrize(dftb.structure)
                    self.structure = dftb.structure
                    dftb.get_geometry()
                    dftb.roll_outputs(irun)
                    dftb.options['CalculateForces'] = True
                    dftb.driver = {}
                    dftb.set_inputs()
                    dftb.run()
                    if self.waiting:
                        dftb.runner.wait()
                    while dftb.runner.poll() is None:
                        dftb.run_status()
                        time.sleep(10)
                    pcm_log.debug('I am alive!!!')
                    filename = dftb.workdir + os.sep + 'detailed.out'
                    forces, stress, total_energy = read_detailed_out(filename=filename)
                    if stress is None or forces is None or total_energy is None:
                        pcm_log.debug('Avoiding a static run relaxing and exiting')
                        dftb.basic_input()
                        dftb.driver['LatticeOpt'] = False
                        # Decreasing the target_forces to avoid the final static
                        # calculation of raising too much the forces after symmetrization
                        dftb.driver['MaxForceComponent'] = 0.9 * self.target_forces
                        dftb.driver['ConvergentForcesOnly'] = False
                        dftb.driver['MaxSteps'] = 10
                        dftb.hamiltonian['MaxSCCIterations'] = 50
                        print dftb.driver
                        dftb.set_inputs()
                        dftb.run()
                        if self.waiting:
                            dftb.runner.wait()
                        while dftb.runner.poll() is None:
                            time.sleep(10)
                        print 'FINAL:', read_detailed_out(filename=filename)
                    break
            else:
                pcm_log.debug('ID: %s' % os.path.basename(self.workdir))
                filename = dftb.workdir + os.sep + 'dftb_stdout.log'
                if os.path.exists(filename):
                    read_dftb_stdout(filename=filename)
                time.sleep(10)
                # else:
                # log.debug('Files: ' + str(os.listdir(dftb.workdir)))
                #    break

    def get_forces_stress_energy(self):

        filename = self.workdir + os.sep + 'OUTCAR'
        if os.path.isfile(filename):
            vo=VaspOutput(filename)
            vo.outcar_parser()
            forces = vo.forces[-1]
            stress = vo.stress[-1]
            total_energy = vo.energy[-1]
        else:
            forces = None
            stress = None
            total_energy = None
        return forces, stress, total_energy
