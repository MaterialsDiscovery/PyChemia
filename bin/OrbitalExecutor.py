#!/usr/bin/env python

import os
import subprocess
import numpy as np
import pychemia

# OPTIONS
USEDMATPU = 25
NSTEP = 50
TOLVRS = 1E-14
TARGET_NRES2 = 1E-12


if __name__ == "__main__":

    # Creating symbolic links to PBS script and abinit.files

    if os.path.lexists('batch.pbs'):
        os.remove('batch.pbs')
    os.symlink('../batch.pbs', 'batch.pbs')

    if os.path.lexists('abinit.files'):
        os.remove('abinit.files')
    os.symlink('../abinit.files', 'abinit.files')

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
        subprocess.call("mpirun -np 6 abinit < abinit.files > abinit.log 2> /dev/null", shell=True)

        # If everything works fine with ABINIT we have abinit.out
        if os.path.isfile('abinit.out'):

            # The final density matrix is build from the output
            dmatpawu=pychemia.population.orbitaldftu.get_final_dmatpawu('abinit.out')
            odmatpawu=np.array(dmatpawu).reshape(-1,5,5)
            
            # Updating dmatpawu from the output back to input
            abiinput['dmatpawu'] = list(odmatpawu.flatten())
            if os.path.isfile('abinit-i_WFK'):
                abiinput['irdwfk'] = 1
            abiinput.write('abinit.in')

        # Renaming logs and setting WFK back to input
        if os.path.isfile('abinit.log'):
            os.rename('abinit.log', 'abinit_%d.log' % i)
        if os.path.isfile('abinit.in'):
            os.rename('abinit.in', 'abinit_%d.in' % i)
        if os.path.isfile('abinit-o_WFK'):
            os.rename('abinit-o_WFK','abinit-i_WFK')

        # Checking if you should accept the current residual
        # Renaming abinit.out
        nres2=1
        if os.path.isfile('abinit.out'):
            abo=pychemia.code.abinit.AbinitOutput('abinit.out')
            nres2=abo.get_energetics()['nres2'][-1]
            os.rename('abinit.out', 'abinit_%d.out' % i)

        if nres2 < TARGET_NRES2:
            break
    wf=open('COMPLETE','w')
    wf.write("%d\n" % i)
    wf.close()
