#!/usr/bin/env python

import pychemia
import numpy as np
import os
import json
import sys
import getopt


def usage(name):
    print("""
NAME

    %s

DESCRIPTION
    Collects information about all the structures found during a Minimal Hoping Method using the original
    code from Amsler,

OPTIONS

    --help, -h
        Return information on the options and use of this script

    --directory, -d <string>
        Directory where the structure search stored the data

    --left, -l <string>
        Specie that will appear on the left side of the convex hull, its composition ratio will be 0.0

    --right, -r <string>
        Specie that will appear on the right side of the convex hull, its composition ratio will be 1.0
""" % os.path.basename(name))


def main(argv):
    try:
        opts, args = getopt.getopt(argv[1:], "hd:l:r:", ["help", "directory=", "left=", "right="])
    except getopt.GetoptError:
        usage(argv[0])
        sys.exit(2)

    if len(opts) == 0:
        usage(argv[0])
        sys.exit(2)

    # Default Values
    directory = None
    left = None
    right = None

    for opt, arg in opts:
        if opt in ("-h", "--help"):
            usage(argv[0])
            sys.exit()
        elif opt in ("-d", "--directory"):
            directory = arg
        elif opt in ("-l", "--left"):
            left = arg
        elif opt in ("-r", "--right"):
            right = arg

    if directory is None or left is None or right is None:
        usage(argv[0])
        sys.exit(1)

    if not os.path.isdir(directory):
        print('Not such directory ', directory)
        sys.exit(1)

    basedir = directory
    specie1 = right
    specie2 = left

    dirs1 = sorted([x for x in os.listdir(basedir) if os.path.isdir(basedir + os.sep + x)])
    ret = []
    for idir in dirs1:
        dirs2 = sorted([x for x in os.listdir(basedir + os.sep + idir) if x[:3] == 'xxx'])
        for idir2 in dirs2:
            path = basedir + os.sep + idir + os.sep + idir2

            if not os.path.exists(path + os.sep + 'CONTCAR') or os.path.getsize(path + os.sep + 'CONTCAR') == 0:
                print("Missing CONTCAR on %s/%s" % (idir, idir2))

                if not os.path.exists(path + os.sep + 'POSCAR') or os.path.getsize(path + os.sep + 'POSCAR') == 0:
                    print("Missing POSCAR on %s/%s" % (idir, idir2))
                    continue
                st = pychemia.code.vasp.read_poscar(path + os.sep + 'POSCAR')
            else:
                st = pychemia.code.vasp.read_poscar(path + os.sep + 'CONTCAR')

            formula = st.formula

            symmetry = pychemia.crystal.CrystalSymmetry(st)
            space_group = symmetry.number(symprec=1e-1)

            vo = pychemia.code.vasp.VaspOutput(path + os.sep + 'OUTCAR')
            if not vo.is_finished:
                continue
            energy = vo.last_energy

            # ana=pychemia.analysis.StructureAnalysis(st)
            # b, c, r = ana.bonds_coordination()

            if vo.forces is not None and len(vo.forces) > 0:
                maxforce = np.max(np.apply_along_axis(np.linalg.norm, 1, vo.forces[-1]))
            else:
                maxforce = 100

            if specie2 not in st.composition:
                ratio = 1
            elif specie1 not in st.composition:
                ratio = 0
            else:
                ratio = float(st.composition[specie1]) / (st.composition[specie2] + st.composition[specie1])

            print(" %30s SPCGRP: %4d   ENERGY_PA: %9.3f   MAXFORCE: %9.2E" % ((idir + os.sep + idir2).ljust(30),
                                                                              space_group, energy / st.natom, maxforce))

            ret.append({'formula': formula,
                        'spcgrp': space_group,
                        'energy': energy,
                        'energy_pa': energy / st.natom,
                        'ratio': ratio,
                        'idir': idir,
                        'idir2': idir2,
                        'natom': st.natom,
                        'maxforce': maxforce,
                        'density': st.density})
    # 'bonds/volume': np.sum(c) / st.volume,
    #                    'max_coordination': np.max(c),
    #                    'coordination': c})

    sorter = []
    for i in ret:
        sorter.append(1E5 * i['ratio'] + i['energy_pa'])

    npsorter = np.array(sorter)
    srt = npsorter.argsort()

    newdata = []
    for i in srt:
        newdata.append(ret[i])

    wf = open('results.json', 'w')
    json.dump(newdata, wf, sort_keys=True, indent=4, separators=(',', ': '))
    wf.close()


if __name__ == "__main__":
    main(sys.argv)
