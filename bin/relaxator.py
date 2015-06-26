#!/usr/bin/env python
import logging
logging.basicConfig(level=logging.DEBUG)
import os
import sys
from pychemia import Structure
from pychemia.code.vasp import ConvergenceKPointGrid, VaspRelaxator, read_poscar, ConvergenceCutOffEnergy
from pychemia.runner import report_cover

try:
    from pychemia.symm import symmetrize
except ImportError:
    def symmetrize(x):
        return x

version = 0.1


def help_info():
    print """
Vasp relaxator with convergence in K-Point Grid and Cut-Off Energy (ENCUT)'
Use:
    relaxator.py [--binary vasp] [ --nparal 8] [--structure_file POSCAR, filename.cif, structure.json]
"""


def cleaner():
    for ifile in ['DOSCAR', 'EIGENVAL', 'IBZKPT', 'OSZICAR', 'PCDAT', 'PROCAR', 'WAVECAR', 'XDATCAR']:
        if os.path.exists(ifile):
            os.remove(ifile)

binary = 'vasp'
workdir = '.'
nparal = 8
structure_file = None
target_forces = 1E-4
energy_tolerance = 5E-3

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
        elif option == 'energy_tolerance':
            energy_tolerance = float(sys.argv[i + 1])
        elif option == 'target_forces':
            target_forces = float(sys.argv[i + 1])
        elif option == 'structure_file':
            structure_file = sys.argv[i + 1]

if structure_file is None:
    help_info()
    exit(1)

st = None
if not os.path.isfile(structure_file):
    print 'Could not find %s' % structure_file
    exit(1)
if structure_file[-4:].lower() == 'json':
    st = Structure.load_json(structure_file)
elif structure_file[-3:].lower() == 'cif':
    import pychemia.external.pymatgen
    st = pychemia.external.pymatgen.cif2structure(structure_file)[0]
elif structure_file[-6:].lower() == 'poscar':
    st = read_poscar(structure_file)
else:
    print "Filename to relax not recognized as ('.cif', '.json' or POSCAR) %s" % structure_file
    exit(1)

print report_cover()

# First Round (Relaxing the original structure)
print '\nFirst Round'
print '===========\n'

cleaner()
print '\nConvergence of K-Point Grid'
print '---------------------------\n'
ck = ConvergenceKPointGrid(st, encut=1.3, nparal=nparal, energy_tolerance=energy_tolerance)
ck.run()
ck.save()
ck.plot()
kp = ck.best_kpoints

cleaner()
print '\nConvergence of Cut-off Energy'
print '-----------------------------\n'
ce = ConvergenceCutOffEnergy(st, energy_tolerance=energy_tolerance, nparal=nparal, kpoints=kp)
ce.run()
ce.save()
ce.plot()
encut = ce.best_encut

os.rename('convergence_encut.json', 'convergence_encut_phase1.json')
os.rename('convergence_encut.pdf', 'convergence_encut_phase1.pdf')
os.rename('convergence_kpoints.json', 'convergence_kpoints_phase1.json')
os.rename('convergence_kpoints.pdf', 'convergence_kpoints_phase1.pdf')

cleaner()
print '\nIonic Relaxation'
print '----------------\n'
vr = VaspRelaxator(workdir=workdir, structure=st,
                   relaxator_params={'encut': encut, 'kp_grid': kp.grid, 'nparal': nparal},
                   target_forces=10*target_forces)
vr.run()

structure = vr.get_final_geometry()
structure.save_json(workdir+os.sep+'structure_phase1.json')

# Second Round (Symetrize structure and redo convergences)
st = symmetrize(structure)

print '\nSecond Round'
print '============'

cleaner()
print '\nConvergence of K-Point Grid'
print '---------------------------\n'
ck = ConvergenceKPointGrid(st, encut=encut, nparal=nparal, energy_tolerance=energy_tolerance)
ck.run()
ck.save()
ck.plot()
kp = ck.best_kpoints

cleaner()
print '\nConvergence of Cut-off Energy'
print '-----------------------------\n'
ce = ConvergenceCutOffEnergy(st, energy_tolerance=energy_tolerance, nparal=nparal, kpoints=kp)
ce.run()
ce.save()
ce.plot()
encut = ce.best_encut

os.rename('convergence_encut.json', 'convergence_encut_phase2.json')
os.rename('convergence_encut.pdf', 'convergence_encut_phase2.pdf')
os.rename('convergence_kpoints.json', 'convergence_kpoints_phase2.json')
os.rename('convergence_kpoints.pdf', 'convergence_kpoints_phase2.pdf')

cleaner()
print '\nIonic Relaxation'
print '----------------\n'
vr = VaspRelaxator(workdir=workdir, structure=st,
                   relaxator_params={'encut': encut, 'kp_grid': kp.grid, 'nparal': nparal},
                   target_forces=target_forces)
vr.run()

structure = vr.get_final_geometry()
structure.save_json(workdir+os.sep+'structure_phase2.json')

cleaner()

os.remove('to_relax')
