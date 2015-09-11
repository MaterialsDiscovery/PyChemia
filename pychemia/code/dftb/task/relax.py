import os
import time
from .._dftb import DFTBplus, read_detailed_out, read_dftb_stdout, read_geometry_gen
from pychemia import pcm_log
from pychemia.dft import KPoints
from pychemia.code import Relaxator

try:
    from pychemia.symm import symmetrize
except ImportError:
    def symmetrize(structure):
        return structure


class Relaxation(Relaxator):

    def __init__(self, structure, relaxator_params=None, workdir='.', kpoints=None, target_forces=1E-3, waiting=False):

        self.workdir = workdir
        self.initial_structure = structure
        self.slater_path = None
        self.symmetrize = False
        if relaxator_params is None:
            relaxator_params = {'slater_path': '.'}

        self.set_params(relaxator_params)

        if self.symmetrize:
            self.initial_structure = symmetrize(structure)
        self.structure = self.initial_structure.copy()
        self.waiting = waiting
        if kpoints is None:
            self.kpoints = KPoints()
            self.kpoints.set_optimized_grid(self.structure.lattice, density_of_kpoints=10000)
        else:
            self.kpoints = kpoints

        Relaxator.__init__(self, target_forces)

    def set_params(self, params):
        assert(isinstance(params, dict))
        if 'slater_path' not in params:
            params['slater_path'] = '.'
        if isinstance(params['slater_path'], basestring):
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
        score = 0
        dftb = DFTBplus()
        dftb.initialize(workdir=self.workdir, structure=self.structure, kpoints=self.kpoints)
        dftb.set_slater_koster(search_paths=self.slater_path)
        dftb.basic_input()
        dftb.driver['LatticeOpt'] = False
        # Decreasing the target_forces to avoid the final static
        # calculation of raising too much the forces after symmetrization
        dftb.driver['MaxForceComponent'] = 0.5 * self.target_forces
        dftb.driver['ConvergentForcesOnly'] = True
        dftb.driver['MaxSteps'] = 50
        dftb.hamiltonian['MaxSCCIterations'] = 50
        dftb.set_inputs()
        print 'Launching DFTB+ with target force of %9.2E ' % dftb.driver['MaxForceComponent']
        dftb.run()
        if self.waiting:
            dftb.runner.wait()
        while True:
            if dftb.runner is not None and dftb.runner.poll() is not None:
                pcm_log.info('Execution completed. Return code %d' % dftb.runner.returncode)
                booleans, geom_optimization, stats=read_dftb_stdout()
                print 'Converged: ', stats['ion_convergence']
                print 'SCC:', geom_optimization['nscc_per_ionstep']
                print 'Forces:', geom_optimization['max_force'][-1]
                if 'max_lattice_force' in geom_optimization:
                    print 'Stress:', geom_optimization['max_lattice_force'][-1]

                filename = dftb.workdir + os.sep + 'detailed.out'
                if not os.path.exists(filename):
                    pcm_log.error('Could not find ' + filename)
                    break
                good_forces, good_stress, good_data = self.relaxation_status()

                if not good_data:
                    # This happens when all the SCC are completed without convergence
                    dftb.driver['ConvergentForcesOnly'] = False
                else:
                    dftb.driver['ConvergentForcesOnly'] = True

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
                    dftb.structure.save_json(dftb.workdir + os.sep + 'structure_current.json')
                    if self.symmetrize:
                        dftb.structure = symmetrize(dftb.structure)
                    self.structure = dftb.structure

                    dftb.get_geometry()
                    dftb.roll_outputs(irun)
                    dftb.set_inputs()
                    irun += 1
                    print 'Launching DFTB+ with target force of %9.2E ' % dftb.driver['MaxForceComponent']
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
                    print 'Launching DFTB+ with static evaluation of forces '
                    dftb.run()
                    if self.waiting:
                        dftb.runner.wait()
                    while dftb.runner.poll() is None:
                        dftb.run_status()
                        time.sleep(10)
                    pcm_log.debug('I am alive!!!')
                    filename = dftb.workdir + os.sep + 'detailed.out'

                    ret = read_detailed_out(filename=filename)
                    forces = ret['forces']
                    stress = ret['stress']
                    total_energy = ret['total_energy']

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

    def quality(self, dftb, score):

        booleans, geom_optimization, stats = read_dftb_stdout(filename=dftb.workdir + os.sep + 'dftb_stdout.log')

        # Increase the score on each iteration
        score += 1

        if stats['ion_convergence']:
            score += 0

        if 'max_force' in geom_optimization:
            max_force = geom_optimization['max_force'][-1]
        else:
            return score

        if 'max_lattice_force' in geom_optimization:
            max_lattice_force = geom_optimization['max_lattice_force'][-1]
        else:
            max_lattice_force = None

        if max_lattice_force is not None and max_force < self.target_forces and max_lattice_force < self.target_forces:
            pcm_log.debug('Target forces and stress achieved (score +100)')
            # Increase for a high value to exit
            score += 100

        if 'max_force' in geom_optimization and geom_optimization['max_force'][-1] > self.target_forces:
            # log.debug('Target forces not achieved (score +1)')
            score += 0

        if 'max_force' in geom_optimization and geom_optimization['max_force'][-1] < geom_optimization['max_force'][0]:
            # log.debug('Forces are decreasing (score +1)')
            score += 0

        # Only for Lattice Optimization
        if 'max_lattice_force' in geom_optimization:
            if geom_optimization['max_lattice_force'][-1] > self.target_forces:
                # log.debug('Target stress not achieved (score +1)')
                score += 0

            if geom_optimization['max_lattice_force'][-1] < geom_optimization['max_lattice_force'][0]:
                # log.debug('Stress is decreasing (score +1)')
                score += 0

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
