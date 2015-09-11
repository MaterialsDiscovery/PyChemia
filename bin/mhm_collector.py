#!/usr/bin/env python

import pychemia
import numpy as np
import os
import json
import sys

basedir = sys.argv[1]
specie1 = sys.argv[2]
specie2 = sys.argv[3]

dirs1 = [x for x in os.listdir(basedir) if os.path.isdir(basedir+os.sep+x)]

ret = []
for idir in dirs1:
    dirs2 = [x for x in os.listdir(basedir+os.sep+idir) if x[:3] == 'xxx']
    for idir2 in dirs2:
        path = basedir+os.sep+idir+os.sep+idir2

        if not os.path.exists(path+os.sep+'CONTCAR') or os.path.getsize(path+os.sep+'CONTCAR') == 0:
            print "Missing CONTCAR on %s/%s" % (idir, idir2)

            if not os.path.exists(path+os.sep+'POSCAR') or os.path.getsize(path+os.sep+'POSCAR') == 0:
                print "Missing POSCAR on %s/%s" % (idir, idir2)	
                continue
            st = pychemia.code.vasp.read_poscar(path+os.sep+'POSCAR')
        else:
            st = pychemia.code.vasp.read_poscar(path+os.sep+'CONTCAR')

        formula = st.formula

        symmetry = pychemia.symm.StructureSymmetry(st)
        space_group = symmetry.number(symprec=1e-1)

        vo = pychemia.code.vasp.VaspOutput(path+os.sep+'OUTCAR')
        if not vo.is_finished:
            continue
        energy = vo.last_energy

        if vo.forces is not None and len(vo.forces) > 0:
            maxforce = np.max(np.apply_along_axis(np.linalg.norm, 1, vo.forces[-1]))
        else:
            maxforce = 100

        if specie2 not in st.composition:
            ratio = 1
        elif specie1 not in st.composition:
            ratio = 0
        else:
            ratio = float(st.composition[specie1])/(st.composition[specie2]+st.composition[specie1])
    
        print "%20s %20s %20s %4d %12.3f %12.3E" % (idir, idir2, formula, space_group, energy, maxforce)

        ret.append({'formula': formula,
                    'spcgrp': space_group,
                    'energy': energy,
                    'energy_pa': energy/st.natom,
                    'ratio': ratio,
                    'idir': idir,
                    'idir2': idir2,
                    'natom': st.natom,
                    'maxforce': maxforce})

sorter = []
for i in ret:
    sorter.append(1E5*i['ratio']+i['energy_pa'])

npsorter = np.array(sorter)
srt = npsorter.argsort()

newdata = []
for i in srt:
    newdata.append(ret[i])

wf = open('results.json', 'w')
json.dump(newdata, wf, sort_keys=True, indent=4, separators=(',', ': '))
wf.close()
