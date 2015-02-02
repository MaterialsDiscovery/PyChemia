__author__ = 'Guillermo Avendano Franco'

import os
import time
import numpy as np
from _dftb import DFTBplus, read_detailed_out, read_dftb_stdout, read_geometry_gen
from pychemia.dft import KPoints
import logging

logging.basicConfig(level=logging.DEBUG)

class Relaxator():

    def __init__(self, workdir, structure, slater_path, target_forces=1E-3):

        self.workdir = workdir
        self.structure = structure
        self.slater_path = slater_path
        self.target_forces = target_forces

    def run(self):

        irun = 0
        score = 10
        dftb = DFTBplus()
        kpoints = KPoints(kmode='gamma', grid=[7, 7, 7])
        dftb.initialize(workdir=self.workdir, structure=self.structure, kpoints=kpoints)
        dftb.set_slater_koster(search_paths=self.slater_path)
        dftb.basic_input()
        dftb.driver['LatticeOpt'] = False
        dftb.driver['MaxForceComponent'] = self.target_forces
        dftb.driver['ConvergentForcesOnly'] = True
        dftb.driver['MaxSteps'] = 50
        dftb.hamiltonian['MaxSCCIterations'] = 50
        dftb.set_inputs()
        dftb.run()
        while True:
            if dftb.runner is not None and dftb.runner.poll() is not None:
                dftb.runner = None
                logging.info('Execution completed!')
                filename = dftb.workdir+os.sep+'detailed.out'
                if not os.path.exists(filename):
                    logging.error('Could not find ' + filename)
                    break
                forces, stress, total_energy = read_detailed_out(filename=filename)

                if forces is None and stress is None and total_energy is None:
                    # This happens when all the SCC are completed without convergence
                    score -= 1
                    dftb.driver['ConvergentForcesOnly'] = False
                else:
                    dftb.driver['ConvergentForcesOnly'] = True

                # Forces are good if we are at least one order of magnitude higher than target
                if forces is not None and np.max(np.abs(forces.flatten())) < 10*self.target_forces:
                    good_forces = True
                else:
                    if forces is not None:
                        pass
                        #logging.info('Forces: '+str(np.max(np.abs(forces.flatten())), 10*self.target_forces))
                    good_forces = False

                # Stress is good if we are at least one order of magnitude higher than target
                if stress is not None and np.max(np.abs(stress.flatten())) < 10*self.target_forces:
                    good_stress = True
                else:
                    if stress is not None:
                        logging.info('Stress: ' + str(np.max(np.abs(stress.flatten())), 10*self.target_forces))
                    good_stress = False

                score = self.quality(dftb, score)
                logging.debug('Present score : '+str(score))
                if score > 0:

                    if good_forces and good_stress:
                        logging.debug('Convergence: Internals + Cell')
                        dftb.driver['MovedAtoms'] = '1:-1'
                        dftb.driver['LatticeOpt'] = True
                    elif not good_forces and good_stress:
                        logging.debug('Convergence: Internals')
                        dftb.driver['LatticeOpt'] = False
                        dftb.driver['MovedAtoms'] = '1:-1'
                    elif good_forces and not good_stress:
                        logging.debug('Convergence: Internals + Cell')
                        dftb.driver['LatticeOpt'] = True
                        dftb.driver['MovedAtoms'] = '1:-1'

                    dftb.structure = read_geometry_gen(dftb.workdir+os.sep+'geo_end.gen')
                    dftb.get_geometry()
                    dftb.roll_outputs(irun)
                    dftb.set_inputs()
                    irun += 1
                    dftb.run()
                else:
                    logging.debug('Not running anymore')
                    break
            else:
                filename = dftb.workdir+os.sep+'dftb_stdout.log'
                if os.path.exists(filename):
                    read_dftb_stdout(filename=filename)
                else:
                    logging.debug('Files: ' + str(os.listdir(dftb.workdir)))
                    break
                time.sleep(30)

    def quality(self, dftb, score):

        booleans, geom_optimization, stats = read_dftb_stdout(filename=dftb.workdir+os.sep+'dftb_stdout.log')

        if stats['ion_convergence']:
            score -= 1

        if 'max_force' in geom_optimization:
            max_force = geom_optimization['max_force'][-1]
        else:
            return score

        if 'max_lattice_force' in geom_optimization:
            max_lattice_force = geom_optimization['max_lattice_force'][-1]
        else:
            max_lattice_force = None

        if max_lattice_force is not None and max_force < self.target_forces and max_lattice_force < self.target_forces:
            logging.debug('Target forces and stress achieved (score -100)')
            score -= 100
        else:
            pass
            #logging.debug('Target Forces: ', self.target_forces)

        if 'max_force' in geom_optimization and geom_optimization['max_force'][-1] > self.target_forces:
            logging.debug('Target forces not achieved (score +1)')
            score += 1

        if 'max_force' in geom_optimization and geom_optimization['max_force'][-1] < geom_optimization['max_force'][0]:
            logging.debug('Forces are decreasing (score +1)')
            score += 1

        # Only for Lattice Optimization
        if 'max_lattice_force' in geom_optimization:
            if geom_optimization['max_lattice_force'][-1] > self.target_forces:
                logging.debug('Target stress not achieved (score +1)')
                score += 1

            if geom_optimization['max_lattice_force'][-1] < geom_optimization['max_lattice_force'][0]:
                logging.debug('Stress is decreasing (score +1)')
                score += 1

        return score
