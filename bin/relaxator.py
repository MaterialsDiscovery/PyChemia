#!/usr/bin/env python

import os
import sys
import logging
logging.basicConfig(level=logging.DEBUG)
import pychemia
from pychemia.code.vasp import Convergence_Kpoints, VaspRelaxator, VaspJob


def help_info():
    print 'Vasp relaxator'

workdir = '.'

if os.path.isfile(workdir+os.sep+'structure1.json'):
    structure = pychemia.Structure.load_json(workdir+os.sep+'structure1.json')
else:
    structure = pychemia.Structure.load_json(workdir+os.sep+'structure.json')

binary = 'vasp'
nparal = 2
version = 0.1

for i in range(1, len(sys.argv)):
    if sys.argv[i].startswith('--'):
        option = sys.argv[i][2:]
        # fetch sys.argv[1] but without the first two characters
        if option == 'version':
            print(version)
            sys.exit()
        elif option == 'help':
            help_info()
            sys.exit()
        elif option == 'binary':
            binary = sys.argv[i + 1]
        elif option == 'nparal':
            nparal = int(sys.argv[i + 1])

# FIRST ROUND

ck = Convergence_Kpoints(workdir, structure, binary=binary, nparal=nparal, energy_tolerance=1E-2, recover=False)
ck.run()

grid = ck.best_grid

kp = ck.best_kpoints()

vj = VaspJob()
vj.initialize(workdir, structure, kp, binary=binary)

rel = VaspRelaxator(workdir, structure, {'kpoints_grid': grid, 'nparal': nparal}, binary=binary)
rel.run()

structure = rel.get_final_geometry()
structure.save_json(workdir+os.sep+'structure1.json')

# SECOND ROUND

ck = Convergence_Kpoints(workdir, structure, binary=binary, nparal=nparal, energy_tolerance=1E-3, recover=True)
ck.run()

grid = ck.best_grid

kp = ck.best_kpoints()

vj = VaspJob()
vj.initialize(workdir, structure, kp, binary=binary)

rel = VaspRelaxator(workdir, structure, {'kpoints_grid': grid, 'nparal': nparal}, binary=binary)
rel.run()

structure = rel.get_final_geometry()
structure.save_json(workdir+os.sep+'structure2.json')
