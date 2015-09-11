#!/usr/bin/env python

import sys
import json
import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull


def spcgrp_props(space_group):
    cst_marker = '.'
    cst_color = '#000000'
    cst_marker_size = 2
    cst_zorder = 2
    cst_fillstyle = 'full'
    if space_group < 16:
        cst_label = 'Triclinic'
    elif space_group < 75:
        cst_marker = 'o'
        cst_label = 'Orthorhombic'
        cst_marker_size = 4
        cst_zorder = 9
        cst_color = '#0000ff'
        cst_fillstyle = 'none'
    elif space_group < 143:
        cst_label = 'Tetragonal'
    elif space_group < 168:
        cst_label = 'Trigonal'
    elif space_group < 195:
        cst_marker = '^'
        cst_label = 'Hexagonal'
        cst_marker_size = 5
        cst_zorder = 8
        cst_color = '#ff00ff'
        cst_fillstyle = 'none'
    else:
        cst_marker = '*'
        cst_label = 'Cubic'
        cst_marker_size = 10
        cst_zorder = 10
        cst_color = '#ff0000'
    return cst_marker, cst_color, cst_label, cst_marker_size, cst_zorder, cst_fillstyle

data = json.load(open('results.json'))
fig = plt.figure(figsize=(11, 8.5))

EA = None
EB = None

for idata in data:
    if idata['ratio'] == 0:
        EA = idata['energy_pa']
        print 'EA', EA
    elif idata['ratio'] == 1:
        EB = idata['energy_pa']
        print 'EB', EB

if EA is None or EB is None:
    print 'Pure elements not found, formation energy cannot be computed'
    sys.exit(1)

points = []
for idata in data:
    x = idata['ratio']
    spcgrp = idata['spcgrp']
    alpha = 0.3
    marker, color, lab, m, z, fs = spcgrp_props(spcgrp)
    y = idata['energy_pa'] - (1-x)*EA - x*EB

    plt.plot(x, y, marker=marker, ms=m, color=color, fillstyle=fs, zorder=z)
    points.append([x, idata['energy_pa'] - (1-x)*EA - x*EB])

points = np.array(points)
hull = ConvexHull(points)

for simplex in hull.simplices:
    if points[simplex, 1][0] <= 0.0 and points[simplex, 1][1] <= 0.0:
        if points[simplex, 1][0] == 0.0 and points[simplex, 1][1] == 0.0:
            continue
        plt.plot(points[simplex, 0], points[simplex, 1], 'k-', zorder=1)
        print points[simplex, 1]

ylims = plt.ylim()

for spcgrp in [74, 194, 230]:
    marker, color, lab, m, z, fs = spcgrp_props(spcgrp)
    print marker
    plt.plot(-100, -100, marker, ms=m, fillstyle=fs, color=color, label=lab)

plt.xlim(-0.05, 1.05)
plt.legend(loc=9, prop={'size': 10}, numpoints=1)
plt.subplots_adjust(left=0.12, bottom=0.13, right=0.98, top=0.96, wspace=None, hspace=None)
plt.xlabel(r'Composition balance')
plt.ylabel(r'Formation Energy [eV]')
plt.savefig('figure.pdf')
