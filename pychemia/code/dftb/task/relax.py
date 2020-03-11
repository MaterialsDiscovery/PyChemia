import os
import time

import numpy as np

from pychemia import pcm_log
from pychemia.code import Relaxator
from pychemia.crystal import KPoints
from ..dftb import DFTBplus, read_detailed_out, read_dftb_stdout, read_geometry_gen

try:
    from pychemia.symm import symmetrize
except ImportError:
    def symmetrize(structure):
        return structure

INITIAL_SCORE = -25


class Relaxation(Relaxator):
    def __init__(self, structure, relaxator_params=None, workdir='.', kpoints=None, target_forces=1E-3, waiting=False,
                 kp_density=10000, forced=True):

        self.workdir = workdir
        self.initial_structure = structure
        self.slater_path = None
        self.symmetrize = False
        self.forced = forced
        if relaxator_params is None:
            relaxator_params = {'slater_path': '.'}

        self.set_params(relaxator_params)

        if self.symmetrize:
            self.initial_structure = symmetrize(structure)
        self.structure = self.initial_structure.copy()
        self.waiting = waiting
        if kpoints is None:
            self.kpoints = KPoints.optimized_grid(self.structure.lattice, kp_density=kp_density)
        else:
            self.kpoints = kpoints

        Relaxator.__init__(self, target_forces)

    def set_params(self, params):
        assert (isinstance(params, dict))
        if 'slater_path' not in params:
            params['slater_path'] = '.'
        if isinstance(params['slater_path'], str):
            assert os.path.exists(params['slater_path'])
            self.slater_path = [params['slater_path']]
        else:
            self.slater_path = params['slater_path']
            try:
                for x in self.slater_path:
                    assert os.path.exists(x)
            except TypeError:
                raise ValueError('Missing a valid slater_path or list of slater_paths')
        if 'symmetrize' in params and params['symmetrize'] is True:
            self.symmetrize = True

    def run(self):

        irun = 0
        score = INITIAL_SCORE
        dftb = DFTBplus(workdir=self.workdir)
        dftb.initialize(structure=self.structure, kpoints=self.kpoints)
        dftb.set_slater_koster(search_paths=self.slater_path)
        dftb.basic_input()
        dftb.driver['LatticeOpt'] = False
        # Decreasing the target_forces to avoid the final static
        # calculation of raising too much the forces after symmetrization
        dftb.driver['MaxForceComponent'] = self.target_forces
        dftb.driver['ConvergentForcesOnly'] = True
        dftb.driver['MaxSteps'] = 100
        dftb.hamiltonian['MaxSCCIterations'] = 20
        dftb.set_inputs()
        print('Launching DFTB+ with target force of %9.2E ' % dftb.driver['MaxForceComponent'])
        dftb.run()
        if self.waiting:
            dftb.runner.wait()
        while True:
            if dftb.runner is not None and dftb.runner.poll() is not None:
                pcm_log.info('Execution completed. Return code %d' % dftb.runner.returncode)
                stdo = read_dftb_stdout(filename=self.workdir + os.sep + 'dftb_stdout.log')

                good_forces, good_stress = self.relaxation_status()
                if 'max_force' in stdo:
                    print('Converged: %s\t Max Force: %9.3e\t MaxForceComponent: %9.3e' % (stdo['ion_convergence'],
                                                                                           stdo['max_force'],
                                                                                           self.target_forces))

                filename = dftb.workdir + os.sep + 'detailed.out'
                if not os.path.exists(filename):
                    pcm_log.error('Could not find ' + filename)
                    break

                if not good_forces and not good_stress:
                    # This happens when all the SCC are completed without convergence
                    dftb.driver['ConvergentForcesOnly'] = False
                else:
                    dftb.driver['ConvergentForcesOnly'] = True

                score = self.quality(score)
                pcm_log.debug('Score :  %d Good Forces: %s   Good Stress: %s' % (score, good_forces, good_stress))
                if score < 0:

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

                    # lets change the positions if the score have lowered to -10
                    if score == -10 and self.forced:
                        dftb.structure.positions += 0.2 * np.random.rand(dftb.structure.natom, 3) - 0.1
                        dftb.structure.positions2reduced()
                        dftb.structure.set_cell(1.1 * dftb.structure.cell)
                    if score == -1 and self.forced:
                        dftb.structure = dftb.structure.random_cell(dftb.structure.composition)
                        print('RANDOM STRUCTURE')
                        print(dftb.structure)
                        score = INITIAL_SCORE

                    dftb.structure.save_json(dftb.workdir + os.sep + 'structure_current.json')
                    if self.symmetrize:
                        dftb.structure = symmetrize(dftb.structure)
                    self.structure = dftb.structure

                    dftb.get_geometry()
                    dftb.roll_outputs(irun)
                    dftb.set_inputs()
                    irun += 1
                    print('Launching DFTB+ with target force of %9.2E ' % dftb.driver['MaxForceComponent'])
                    dftb.run()
                    if self.waiting:
                        dftb.runner.wait()
                else:
                    pcm_log.debug('Final static calculation')
                    dftb.structure = self.get_final_geometry()
                    dftb.structure.save_json(dftb.workdir + os.sep + 'structure_final.json')
                    if self.symmetrize:
                        dftb.structure = symmetrize(dftb.structure)
                    self.structure = dftb.structure
                    dftb.get_geometry()
                    dftb.roll_outputs(irun)
                    dftb.options['CalculateForces'] = True
                    dftb.driver = {}
                    dftb.set_inputs()
                    print('Launching DFTB+ with static evaluation of forces ')
                    dftb.run()
                    if self.waiting:
                        dftb.runner.wait()
                    while dftb.runner.poll() is None:
                        dftb.run_status()
                        time.sleep(10)
                    print('Completed Static run')
                    forces, stress, total_energy = self.get_forces_stress_energy()

                    if stress is None or forces is None or total_energy is None:
                        pcm_log.debug('Null Forces, Stress or Energy, relaxing and exiting')
                        dftb.basic_input()
                        dftb.driver['LatticeOpt'] = False
                        # Decreasing the target_forces to avoid the final static
                        # calculation of raising too much the forces after symmetrization
                        dftb.driver['MaxForceComponent'] = 0.9 * self.target_forces
                        dftb.driver['ConvergentForcesOnly'] = False
                        dftb.driver['MaxSteps'] = 10
                        dftb.hamiltonian['MaxSCCIterations'] = 50
                        print(dftb.driver)
                        dftb.set_inputs()
                        dftb.run()
                        if self.waiting:
                            dftb.runner.wait()
                        while dftb.runner.poll() is None:
                            time.sleep(10)
                        print('FINAL:', read_detailed_out(filename=filename))
                        forces, stress, total_energy = self.get_forces_stress_energy()
                        if stress is None or forces is None or total_energy is None:
                            pcm_log.debug('Again Null Forces, Stress or Energy, Randomizing Structure')
                            dftb.structure = dftb.structure.random_cell(dftb.structure.composition)
                            print('RANDOM STRUCTURE')
                            print(dftb.structure)
                            score = INITIAL_SCORE
                        else:
                            break
                    else:
                        break
            else:
                pcm_log.debug('ID: %s' % os.path.basename(self.workdir))
                filename = dftb.workdir + os.sep + 'dftb_stdout.log'
                if os.path.exists(filename):
                    stdo = read_dftb_stdout(filename=filename)
                    print('Number of steps:', len(stdo['Geometry_Steps']))
                    if len(stdo['Geometry_Steps']) > 1:
                        line = 'Energy behavior: '
                        prev_energy = stdo['Geometry_Steps'][0]['Total Energy']['value']
                        line += ' %7.3f ' % prev_energy
                        for step in stdo['Geometry_Steps'][1:]:
                            new_energy = step['Total Energy']['value']
                            if prev_energy > new_energy:
                                line += '>'
                            else:
                                line += '<'
                            prev_energy = new_energy
                        finene = stdo['Geometry_Steps'][-1]['Total Energy']['value']
                        line += ' %7.3f' % finene
                        print(line)
                time.sleep(10)

    def quality(self, score):

        good_forces, good_stress = self.relaxation_status()
        if good_forces and good_stress:
            print('Finished with forces and stress under target_forces')
            score = 0
        elif good_forces:
            print('Finished with forces under target_forces (not stress)')
            score = score
        else:
            # Increase the score on each iteration
            score += 1

        if self.structure.density < 0.1:
            print('Very small density = Bad Structure')
            score = -1

        return score

    def get_forces_stress_energy(self):

        filename = self.workdir + os.sep + 'detailed.out'
        if os.path.isfile(filename):
            ret = read_detailed_out(filename=filename)
            forces = ret['forces']
            stress = ret['stress']
            total_energy = ret['total_energy']
        else:
            forces = None
            stress = None
            total_energy = None
        return forces, stress, total_energy

    def get_final_geometry(self):
        geometry = self.workdir + os.sep + 'geo_end.gen'
        if os.path.isfile(geometry):
            return read_geometry_gen(geometry)
        else:
            # For static calculations the 'geo_end.gen' is not generated
            # returning the internal structure
            return self.structure
