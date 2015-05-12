__author__ = 'Guillermo Avendano Franco'

import os
import shutil
import pychemia
import time
import numpy as np
from pychemia import pcm_log
from pychemia.dft import KPoints

try:
    from pychemia.symm import symmetrize
except ImportError:
    def symmetrize(structure):
        return structure
from pychemia.code import Relaxator


class RelaxPopulation():

    def __init__(self, population, basedir):
        self.population = population
        self.basedir = basedir
        self.AbinitJobs = {}

    def create_dirs(self, clean=False):
        if not os.path.isdir(self.basedir):
            os.makedirs(self.basedir)
        elif clean:
            for i in os.listdir(self.basedir):
                path = self.basedir+os.sep+i
                if os.path.isfile(path):
                    os.remove(path)
                elif os.path.isdir(path):
                    shutil.rmtree(path)
        for i in self.population.pcdb.entries.find():
            name = self.basedir+os.sep+str(i['_id'])
            if not os.path.isdir(name):
                os.mkdir(name)

    def create_inputs(self):
        kpoints = pychemia.dft.KPoints(kmode='gamma', grid=[4, 4, 4])
        for i in self.population.pcdb.entries.find():
            name = str(i['_id'])
            workdir = self.basedir+os.sep+name
            struct = pychemia.Structure.from_dict(i)
            job = pychemia.code.abinit.AbinitJob(struct, workdir)
            job.set_kpoints(kpoints)
            job.ion_relax(tolmxf=1E-4, tolrff=1E-2,)
            job.energetics()
            job.write_all()
            self.AbinitJobs[name] = job


class AbinitRelaxator(Relaxator):
    def __init__(self, workdir, structure, relaxator_params, target_forces=1E-3, waiting=False):

        self.workdir = workdir
        self.initial_structure = symmetrize(structure)
        self.structure = self.initial_structure.copy()
        self.psp_path = None
        self.set_params(relaxator_params)
        self.waiting = waiting
        Relaxator.__init__(self, target_forces)

    def set_params(self, params):
        assert(isinstance(params, dict))
        if 'psp_path' not in params:
            raise ValueError('A least one path must to be present to search for Pseudo-potential files')
        else:
            if isinstance(params['psp_path'], basestring):
                assert os.path.exists(params['psp_path'])
                self.psp_path = [params['psp_path']]
            else:
                self.psp_path = params['psp_path']
                try:
                    for x in self.psp_path:
                        assert os.path.exists(x)
                except TypeError:
                    raise ValueError('Missing a valid psp_path or list of psp_paths')

    def run(self):

        irun = 0
        score = 0
        job = pychemia.code.abinit.AbinitJob(struct, workdir)
        kpoints = KPoints(kmode='gamma', grid=[7, 7, 7])
        job.initialize(workdir=self.workdir, structure=self.structure, kpoints=kpoints)
        job.set_psp(search_paths=self.slater_path)
        job.basic_input()
        dftb.driver['LatticeOpt'] = False
        # Decreasing the target_forces to avoid the final static
        # calculation of raising too much the forces after symmetrization
        dftb.driver['MaxForceComponent'] = 0.5 * self.target_forces
        dftb.driver['ConvergentForcesOnly'] = True
        dftb.driver['MaxSteps'] = 50
        dftb.hamiltonian['MaxSCCIterations'] = 50
        dftb.set_inputs()
        dftb.run()
        if self.waiting:
            dftb.runner.wait()
        while True:
            if dftb.runner is not None and dftb.runner.poll() is not None:
                pcm_log.info('Execution completed. Return code %d' % dftb.runner.returncode)
                filename = dftb.workdir + os.sep + 'detailed.out'
                if not os.path.exists(filename):
                    pcm_log.error('Could not find ' + filename)
                    break
                forces, stress, total_energy = read_detailed_out(filename=filename)

                if forces is None or stress is None:
                    # This happens when all the SCC are completed without convergence
                    dftb.driver['ConvergentForcesOnly'] = False
                else:
                    dftb.driver['ConvergentForcesOnly'] = True

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
            forces, stress, total_energy = read_detailed_out(filename=filename)
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

