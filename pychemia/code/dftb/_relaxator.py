__author__ = 'Guillermo Avendano Franco'

import os
import time
import numpy as np
from _dftb import DFTBplus, read_detailed_out, read_dftb_stdout, read_geometry_gen
from pychemia.dft import KPoints

class Relaxator():

    def __init__(self, workdir, structure, slater_path, target_forces=1E-3, target_stress=1E-3):

        self.workdir = workdir
        self.structure = structure
        self.slater_path = slater_path
        self.target_forces = target_forces
        self.target_stress = target_stress


    def run(self):

        irun = 0
        dftb = DFTBplus()
        kpoints = KPoints(kmode='gamma', grid=[7,7,7])
        dftb.initialize(workdir=self.workdir, structure=self.structure, kpoints=kpoints)
        dftb.set_slater_koster(search_paths=self.slater_path)
        dftb.basic_input()
        dftb.driver['LatticeOpt'] = False
        dftb.driver['MaxForceComponent'] = self.target_forces
        dftb.driver['ConvergentForcesOnly'] = True
        dftb.set_inputs()
        dftb.run()
        while True:
            if dftb.runner is not None and dftb.runner.poll() is not None:
                dftb.runner = None
                print 'Execution completed!'
                filename = dftb.workdir+os.sep+'detailed.out'
                if not os.path.exists(filename):
                    print 'Could not find ', filename
                    break
                forces, stress, total_energy = read_detailed_out(filename=filename)
                # Forces are good if we are at least one order of magnitude than target
                good_forces = False
                if forces is not None and np.max(np.abs(forces.flatten())) < 10*self.target_forces:
                    print 'Forces:\n'+str(forces)
                    good_forces = True
                else:
                    good_forces = False

                # Stress is good if we are at least one order of magnitude higher than target
                good_stress = False
                if stres is not None and np.max(np.abs(stress.flatten())) < 10*self.target_stress:
                    print 'Stress:\n'+str(stress)
                    good_stress = True
                else:
                    good_stress = False

                print 'Total Energy:\n'+str(total_energy)

                score = self.quality(dftb)
                if score > 0:
                    print 'Executing once more'

                    if good_forces and good_stress:
                        dftb.driver['MovedAtoms'] = '1:-1'
                        dftb.driver['LatticeOpt'] = True
                    elif good_forces and not good_stress:
                        dftb.driver['LatticeOpt'] = False


                    structure = read_geometry_gen(dftb.workdir+os.sep+'geo_end.gen')
                    print structure
                    dftb.structure = structure
                    dftb.get_geometry()
                    dftb.roll_outputs(irun)
                    dftb.set_inputs()
                    irun += 1
                    dftb.run()
                else:
                    print 'Not running anymore'
                    break
            else:
                filename = dftb.workdir+os.sep+'dftb_stdout.log'
                if os.path.exists(filename):
                    print 'Status: '
                    read_dftb_stdout(filename=filename)
                else:
                    print 'Files: ', os.listdir(dftb.workdir)
                    break
                print 'Sleeping 10 seconds...'
                time.sleep(10)


    def quality(self, dftb):

        booleans, geom_optimization, stats = read_dftb_stdout(filename= dftb.workdir+os.sep+'dftb_stdout.log')
        score = 0

        if 'max_force' not in geom_optimization:
            print 'No forces recorded. Returning 1'
            return 1

        if ('max_force' in geom_optimization and geom_optimization['max_force'][-1] < self.target_forces) and \
                ('max_lattice_force' in geom_optimization and ['max_lattice_force'][-1] < self.target_stress):
            print 'Target forces and stress achieved (score -100)'
            score -= 100

        if 'max_force' in geom_optimization and geom_optimization['max_force'][-1] > self.target_forces:
            print 'Target forces not achieved (score +1)'
            score += 1

        if geom_optimization['max_force'][-1] < geom_optimization['max_force'][0]:
            print 'Forces are decreasing (score +1)'
            score += 1

        # Only for Lattice Optimization
        if 'max_lattice_force' in geom_optimization:
            if geom_optimization['max_lattice_force'][-1] > self.target_stress:
                print 'Target stress not achieved (score +1)'
                score += 1

            if geom_optimization['max_lattice_force'][-1] < geom_optimization['max_lattice_force'][0]:
                print 'Stress is decreasing (score +1)'
                score += 1

        print 'The score is: ', score
        return score
