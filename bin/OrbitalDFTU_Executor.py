#!/usr/bin/env python

from __future__ import print_function
import os
import shutil
import subprocess
import numpy as np
import pychemia
import time

# OPTIONS
USEDMATPU = 23
NSTEP = 46
TOLVRS = 1E-14
TARGET_NRES2 = 1E-12
NRUNS=10

if __name__ == "__main__":
    
    print("ABINIT Orbital Evaluator")
    print("========================\n")

    walltime = os.getenv('PBS_WALLTIME')
    if walltime is None:
        print("Env Variable $PBS_WALLTIME is not present, execute this with a PBS script")
        exit(1)
    walltime = int(walltime)

    nodefile=os.getenv('PBS_NODEFILE')
    if nodefile is None:
        print("Env Variable $PBS_NODEFILE is not present, execute this with a PBS script")
        exit(1)

    print("Nodefile: %s" % nodefile)
    print("Walltime: %d seconds = %d minutes = %d hours)" % (walltime,int(walltime/60), int(walltime/3600)))
    rf=open(nodefile)
    nparal = len(rf.readlines())
    print("Number of cores for MPI: %d" % nparal)
    print("Checking abinit.in...", end='')
    if not os.path.isfile('abinit.in'):
        print('No')
        raise ValueError("ERROR: Could not open abinit.in")
    else:
        print('Yes')

    # Getting the current time, use to compute the remaining time in execution
    start_time=time.time()

    abi = pychemia.code.abinit.AbinitInput('abinit.in')
    print("Checking that abinit.in contains value for dmatpawu...", end='')
    if not 'dmatpawu' in abi.variables:
        print('No')
        raise ValueError("ERROR: Could not open abinit.in")
    else:
        print('Yes')

    print("Checking that abinit.in contains value for lpawu...", end='')
    if not 'lpawu' in abi.variables:
        print('No')
        raise ValueError("ERROR: Could not open abinit.in")
    else:
        print('Yes, max lpawu=%d' % max(abi['lpawu']))

    print('Setting variables for usedmatpu, nstep and tolvrs')
    abi['usedmatpu'] = USEDMATPU
    abi['nstep'] = NSTEP
    abi['tolvrs'] = TOLVRS
    print('Writting modified abinit.in')
    abi.write('abinit.in')

    # Getting the index from the last execution and adding one for the next run
    index = 0
    while True:
        if os.path.isfile('abinit_%d.in' % index):
            print("Found abinit_%d.in, moving to next run" % index)
            index+=1
        else:
            break

    print("Executing run with index: %d" % index)
    while True:
        print("\n")
        print('ABINIT execution %d of %d' % (index+1, NRUNS))
        abi = pychemia.code.abinit.AbinitInput('abinit.in')

        # If possible set the WFK from the output back to input 
        if os.path.isfile('abinit-i_WFK'):
            abi['irdwfk'] = 1
            abi.write('abinit.in')

        # Calling ABINIT
        command_line="mpirun -np %d abinit < abinit.files > abinit.log 2> abinit.err" % nparal
        print('Running; %s' % command_line)
        start_run=time.time()
        subprocess.call(command_line, shell=True)
        end_run=time.time()
        runtime=end_run-start_run
        print('Execution finished, execution took %d minutes' % int(runtime/60))

        if os.path.isfile('abinit.in'):
            shutil.copy2('abinit.in', 'abinit_%d.in' % index)

        # If everything works fine with ABINIT we have abinit.out
        if not os.path.isfile('abinit.out'):
            raise ValueError('File not found: abinit.out')

        print("Reading the output abinit.out...")
        # The final density matrix is build from the outputi
        ndim = 2*max(abi['lpawu'])+1
        try:
            newdmatpawu = pychemia.population.orbitaldftu.get_final_dmatpawu('abinit.out')
            print('New dmatpawu found, number of elements: %d' % len(dmatpawu))
        except:
            newdmatpawu = None
            print("Could not get final dmatpawu from abinit.out")
        if newdmatpawu is not None:
            dmatpawu = newdmatpawu
        else:
            dmatpawu = abi['dmatpawu']

        print('Reshaping to %d matrices %d X %d' % (len(dmatpawu)/(ndim*ndim), ndim, ndim))
        odmatpawu = np.array(dmatpawu).reshape(-1, ndim, ndim)
        params=pychemia.population.orbitaldftu.dmatpawu2params(dmatpawu, ndim)
        print("New parameters obtained for %d matrices" % params['num_matrices'])

        # Updating dmatpawu from the output back to input
        abi['dmatpawu'] = list(odmatpawu.flatten())
        if os.path.isfile('abinit-i_WFK'):
            abi['irdwfk'] = 1
        abi.write('abinit.in')

        # Renaming logs and setting WFK back to input
        if os.path.isfile('abinit.log'):
            os.rename('abinit.log', 'abinit_%d.log' % index)
        if os.path.isfile('abinit-o_WFK'):
            os.rename('abinit-o_WFK','abinit-i_WFK')

        # Checking if you should accept the current residual
        # Renaming abinit.out
        nres2 = 1.0
        if os.path.isfile('abinit.out'):
            abo = pychemia.code.abinit.AbinitOutput('abinit.out')
            nres2 = abo.get_energetics()['nres2'][-1]
            os.rename('abinit.out', 'abinit_%d.out' % index)

        if nres2 < TARGET_NRES2:
            break

        # Current time
        curtime=time.time()
        if curtime+runtime > start_time + walltime:
            print("Based on previous run, it is unlikely that next run will have time to complete, exiting")
            break
        else:
            print("Remaining time %d minutes, time for one more run" % int((start_time + WALLTIME - curtime)/60) )

    wf = open('COMPLETE','w')
    wf.write("%d\n" % index)
    wf.close()
