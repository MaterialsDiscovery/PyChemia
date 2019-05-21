#!/usr/bin/env python

import os
import sys
import json
import pychemia

filename = None
if len(sys.argv) == 2:
    filename = sys.argv[1]
else:
    print('Enter the JSON filename to store the data')
    exit(1)

dirs = [x for x in os.listdir('.') if os.path.isdir(x)]

ret = []
for idir in dirs:
    if os.path.isfile(idir + '/POSCAR'):
        try:
            st = pychemia.code.vasp.read_poscar(idir + '/POSCAR')
        except ValueError:
            print('Bad POSCAR\n\n' + open(idir + '/POSCAR').read())
            continue
        # shutil.copy2(idir+'/POSCAR',idir+'_POSCAR')
        print(idir, st.natom)
    else:
        st = pychemia.Structure.load_json(idir + '/structure.json')
        # shutil.copy2(idir+'/structure.json',idir+'_structure.json')
        print('ERROR:', idir, st.natom)
        continue

    if os.path.isfile(idir + '/OUTCAR'):
        try:
            vo = pychemia.code.vasp.VaspOutput(idir + '/OUTCAR')
        except ValueError:
            print('Error reading Vasp Output @ ' + idir + '/OUTCAR')
            continue
        if not vo.has_forces_stress_energy():
            print('Error extracting forces @ ' + idir + '/OUTCAR')
            continue
    else:
        print('No OUTCAR found @ ' + idir)
        continue

    spacegroup = pychemia.crystal.CrystalSymmetry(st).number()
    energy_pa = vo.final_data['energy']['free_energy'] / st.natom

    data = {'id': idir, 'energy_pa': energy_pa, 'natom': st.natom, 'spacegroup': spacegroup,
            'forces': vo.relaxation_info()['avg_force'],
            'stress': max(vo.relaxation_info()['avg_stress_diag'], vo.relaxation_info()['avg_stress_non_diag'])}

    ret.append(data)

wf = open(filename, 'w')
json.dump(ret, wf, sort_keys=True, indent=4, separators=(',', ': '))
wf.close()
