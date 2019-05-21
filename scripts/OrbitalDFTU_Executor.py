#!/usr/bin/env python

import os
import argparse
import subprocess
import time
import shutil
import socket
import numpy as np
from pychemia.code.abinit import AbinitInput, AbinitOutput
from pychemia.population.orbitaldftu import dmatpawu2params


def which(program):
    import os

    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Orbital DFTU Executor, execute several ABINIT runs changing the "
                                                 "value from dmatpawu from the previous output as input until a "
                                                 "tolerance is reached. The executor assumes the existance of a "
                                                 "command 'abinit' and 'mpirun' and files 'abinit.in' and "
                                                 "'abinit.files' the canonical names for those files in ABINIT")

    parser.add_argument('--usedmatpu', type=int, help='ABINIT Variable usedmatpu for each run (default: 25)', 
                        required=False, default=25, metavar='<N>')
    parser.add_argument('--nstep', type=int, help='ABINIT Variable nstep for each run (default: 50)', 
                        required=False, default=50, metavar='<N>')
    parser.add_argument('--tolvrs', type=float, help='ABINIT Variable tolvrs for each run (default: 1E-14)', 
                        required=False, default=1E-14, metavar='<X>')
    parser.add_argument('--target_nres2', type=float, help='Stopping criteria for this executor (default: 1E-12)', 
                        required=False, default=1E-12, metavar='<X>')
    parser.add_argument('--max_nruns', type=int, help='Maximum number of runs allowed (default: 10)', 
                        required=False, default=10, metavar='<N>')
    parser.add_argument('--nhours', type=int, help='Maximun number of hours, ignored if running through a queue '
                                                   'system (PBS), mandatory otherwise',
                        required=False, default=0, metavar='<N>')
    parser.add_argument('--nparal', type=int, help='Number of cores for use with MPI, ignored if running through a '
                                                   'queue system (PBS), mandatory otherwise',
                        required=False, default=0, metavar='<N>')

    args = parser.parse_args()

    print(" ABINIT Orbital DFT+U Executor")
    print(" =============================\n")

    print("Running on: %s" % socket.gethostname())
    print("Current load: %f %f %f" % os.getloadavg())

    # Checking the consistency of all arguments:
    if args.target_nres2 <= args.tolvrs:
        print("Target value must be bigger than ABINIT internal criteria for tolvrs")
        parser.print_help()
        exit(1)
    if args.usedmatpu >= args.nstep:
        print("Total number of SCF steps 'nstep' must be bigger than 'usedmatpu' the number of steps with 'dmatpawu' "
              "fixed")
        parser.print_help()
        exit(1)

    if not os.path.exists('abinit.in'):
        with open('ERROR', 'w') as wf:
            wf.write('No abinit.in')
        raise RuntimeError("File 'abinit.in' could not be found or its symbolic link is broken")
    if not os.path.exists('abinit.files'):
        with open('ERROR', 'w') as wf:
            wf.write('No abinit.files')
        raise RuntimeError("File 'abinit.files' could not be found or its symbolic link is broken")

    # Checking the existance of "mpirun" and "abinit"
    ret = which('mpirun')
    if ret is None:
        with open('ERROR', 'w') as wf:
            wf.write('No mpirun executable')
        raise RuntimeError("Command 'mpirun' could not be found, maybe you need to load the module first")
    print("mpirun: %s" % ret)    
    ret = which('abinit')
    if ret is None:
        with open('ERROR', 'w') as wf:
            wf.write('No abinit executable')
        raise RuntimeError("Command 'abinit' could not be found, maybe you need to load the module first")
    print("abinit: %s" % ret)

    usedmatpu = args.usedmatpu
    nstep = args.nstep
    tolvrs = args.tolvrs
    target_nres2 = args.target_nres2
    max_nruns = args.max_nruns

    nodefile = os.getenv('PBS_NODEFILE')
    nparal = 1
    if nodefile is not None:
        print("Nodefile: %s" % nodefile)
        rf = open(nodefile)
        nparal = len(rf.readlines())
    elif args.nparal > 0:
        nparal = args.nparal
    else:
        print("ERROR: No queue system detected and no positive value for 'nparal'")
        exit(1)

    walltime = os.getenv('PBS_WALLTIME')
    if walltime is not None:
        walltime = int(walltime)
    elif args.nhours > 0:
        walltime = int(args.nhours*3600)
    else:
        print("ERROR: No queue system detected and no positive value for 'nhours'")
        exit(1)

    print("Walltime: %d seconds = %d minutes = %d hours)" % (walltime, int(walltime/60), int(walltime/3600)))
    print("Number of cores for MPI: %d" % nparal)

    # Getting the current time, use to compute the remaining time in execution
    start_time = time.time()

    abi = AbinitInput('abinit.in')
    print("Checking that abinit.in contains value for dmatpawu...", end='')
    if 'dmatpawu' not in abi.variables:
        print('No')
        raise ValueError("ERROR: Could not open abinit.in")
    else:
        print('Yes')

    print("Checking that abinit.in contains value for lpawu...", end='')
    if 'lpawu' not in abi.variables:
        print('No')
        raise ValueError("ERROR: Could not open abinit.in")
    else:
        print('Yes, max lpawu=%d' % max(abi['lpawu']))

    print('Setting ABINIT variables usedmatpu=%d, nstep=%d and tolvrs=%e' % (usedmatpu, nstep, tolvrs))
    abi['usedmatpu'] = usedmatpu
    abi['nstep'] = nstep
    abi['tolvrs'] = tolvrs
    print('Writting modified abinit.in')
    abi.write('abinit.in')

    # Getting the index from the last execution and adding one for the next run
    index = 0
    while True:
        if os.path.isfile('abinit_%02d.in' % index):
            index += 1
        else:
            break
    
    if index >= max_nruns:
        print("Total number of runs has been achieve already, increse 'max_nruns' if you want to continue")
        parser.print_help()
        exit(1)

    print("Executing run with index: %d" % index)
    while index < max_nruns:
        print('\nABINIT execution %d of %d' % (index+1, max_nruns))
        abi = AbinitInput('abinit.in')

        # If possible set the WFK from the output back to input 
        if os.path.isfile('abinit-i_WFK'):
            abi['irdwfk'] = 1
            abi.write('abinit.in')

        # Calling ABINIT
        command_line = "mpirun -np %d abinit < abinit.files > abinit.log 2> abinit.err" % nparal
        print('Running; %s' % command_line)
        start_run = time.time()
        subprocess.call(command_line, shell=True)
        end_run = time.time()

        # Delete the error file if empty
        if os.path.isfile('abinit.err') and os.path.getsize('abinit.err') == 0:
            os.remove('abinit.err')

        runtime = end_run-start_run
        print('Execution finished, execution took %d minutes' % int(runtime/60))

        # If everything works fine with ABINIT we have abinit.out
        # Otherwise is better to stop the entire run
        if not os.path.isfile('abinit.out'):
            with open('ERROR', 'w') as wf:
                wf.write('No abinit.out')
            raise ValueError('File not found: abinit.out')

        if not os.path.isfile('abinit.log'):
            with open('ERROR', 'w') as wf:
                wf.write('No abinit.log')
            raise ValueError('File not found: abinit.log')

        if os.path.isfile('abinit.in'):
            shutil.copy2('abinit.in', 'abinit_%02d.in' % index)

        # Renaming logs and setting WFK back to input
        if os.path.isfile('abinit.log'):
            print("Renaming abinit.log")
            os.rename('abinit.log', 'abinit_%02d.log' % index)
        if os.path.isfile('abinit.err'):
            print("Renaming abinit.err")
            os.rename('abinit.err', 'abinit_%02d.err' % index)
        if os.path.isfile('abinit-o_WFK'):
            print("Renaming abinit-o_WFK")
            os.rename('abinit-o_WFK', 'abinit-i_WFK')
        if os.path.isfile('abinit.out'):
            print("Renaming abinit.out")
            shutil.copy2('abinit.out', 'abinit_%02d.txt' % index)
            os.rename('abinit.out', 'abinit_%02d.out' % index)
            abiout = 'abinit_%02d.out' % index
        else:
            print("Could not find abinit.out")

        # Opening the output file
        print("Reading the output from 'abinit.out'...")
        abo = AbinitOutput(abiout)
        if not abo.is_finished:
            print("abinit.out is truncated, discarting that output redoing the calculation")
            continue

        # The final density matrix is build from the outputi
        ndim = 2*max(abi['lpawu'])+1

        params = dmatpawu2params(abi['dmatpawu'], ndim)
        print("Euler angles from 'abinit.in':")
        for i in params['euler_angles']:
            for j in i:
                print(" %7.3f" % j, end=' ')
            print("\n", end='')

        dmatpawu = abo.get_dmatpawu()
        if dmatpawu is not None:
            params = dmatpawu2params(dmatpawu, ndim)
            print("Euler angles from 'abinit.out':")
            for i in params['euler_angles']:
                for j in i:
                    print(" %7.3f" % j, end=' ')
                print("\n", end='')
        else:
            print("Could not get final dmatpawu from abinit.out")
            dmatpawu = np.array(abi['dmatpawu']).reshape((-1, ndim, ndim))

        print('Shape of dmatpawu: %s' % str(dmatpawu.shape))

        # Updating dmatpawu from the output back to input
        abi['dmatpawu'] = list(dmatpawu.flatten())
        if os.path.isfile('abinit-i_WFK'):
            abi['irdwfk'] = 1
        abi.write('abinit.in')

        # Checking if you should accept the current residual
        energetics = abo.get_energetics()
        if energetics is None:
            raise RuntimeError("Could not get energetics from the output")
        nres2 = 1.0
        nres2 = energetics['nres2'][-1]

        if nres2 < target_nres2 or index == max_nruns-1:
            wf = open('COMPLETE', 'w')
            wf.write("%d\n" % index)
            wf.close()
            break

        # Current time
        curtime = time.time()
        if curtime+runtime > start_time + walltime:
            print("Based on previous run, it is unlikely that next run will have time to complete, exiting")
            print("The walltime for this job is %d minutes and we have been running for %d minutes" %
                  (int(walltime/60), int((curtime-start_time)/60)))
            break
        else:
            print("Remaining time %d minutes, time for one more run" % int((start_time + walltime - curtime)/60))

        # Incresing index for next run
        index += 1
