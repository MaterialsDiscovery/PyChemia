#!/usr/bin/env python

import os
import shutil
import subprocess
import numpy as np
import pychemia

# OPTIONS
USEDMATPU = 25
NSTEP = 50
TOLVRS = 1E-14
TARGET_NRES2 = 1E-12
NPARAL = 6


if __name__ == "__main__":

    # Editing abinit.in to ensure the OPTIONS
    abiinput = pychemia.code.abinit.InputVariables('abinit.in')
    abiinput['usedmatpu'] = USEDMATPU
    abiinput['nstep'] = NSTEP
    abiinput['tolvrs'] = TOLVRS
    abiinput.write('abinit.in')        

    for i in range(10):

        abiinput = pychemia.code.abinit.InputVariables('abinit.in')

        # If possible set the WFK from the output back to input 
        if os.path.isfile('abinit-i_WFK'):
            abiinput['irdwfk'] = 1
            abiinput.write('abinit.in')

        # Calling ABINIT
        subprocess.call("mpirun -np %d abinit < abinit.files > abinit.log 2> abinit.err" % NPARAL, shell=True)

        if os.path.isfile('abinit.in'):
            shutil.copy2('abinit.in', 'abinit_%d.in' % i)

        # If everything works fine with ABINIT we have abinit.out
        if os.path.isfile('abinit.out'):

            # The final density matrix is build from the output
            dmatpawu = pychemia.population.orbitaldftu.get_final_dmatpawu('abinit.out')
            odmatpawu = np.array(dmatpawu).reshape(-1, 5, 5)
            
            # Updating dmatpawu from the output back to input
            abiinput['dmatpawu'] = list(odmatpawu.flatten())
            if os.path.isfile('abinit-i_WFK'):
                abiinput['irdwfk'] = 1
            abiinput.write('abinit.in')

        # Renaming logs and setting WFK back to input
        if os.path.isfile('abinit.log'):
            os.rename('abinit.log', 'abinit_%d.log' % i)
        if os.path.isfile('abinit-o_WFK'):
            os.rename('abinit-o_WFK', 'abinit-i_WFK')

        # Checking if you should accept the current residual
        # Renaming abinit.out
        nres2 = 1.0
        if os.path.isfile('abinit.out'):
            abo = pychemia.code.abinit.AbinitOutput('abinit.out')
            nres2 = abo.get_energetics()['nres2'][-1]
            os.rename('abinit.out', 'abinit_%d.out' % i)

        if nres2 < TARGET_NRES2:
            break
    wf = open('COMPLETE', 'w')
    wf.write("%d\n" % i)
    wf.close()
