#!/usr/bin/env python

import os
import sys
import json
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull
import getopt


def usage(name):
    print("""
NAME

    %s

DESCRIPTION
    Make a Convex Hull from the data stored in the collected 'results.json' file

OPTIONS

    --help, -h
        Return information on the options and use of this script

    --ylim_bottom, -b <float>
        Lower limit on y axis

    --ylim_top, -t <float>
        Higher limit on y axis

    --output_file, -o <string> (Default: ConvexHull.pdf)
        Name of the figure that will be created

""" % os.path.basename(name))


def spcgrp_props(space_group):
    cst_marker = '.'
    cst_color = '#000000'
    cst_marker_size = 2
    cst_zorder = 2
    cst_fillstyle = 'full'
    if space_group < 3:
        cst_label = 'Triclinic'
    elif space_group < 16:
        cst_label = 'Monoclinic'
        cst_marker = '>'
        cst_marker_size = 3
        cst_zorder = 8
        cst_color = '#00aaaa'
    elif space_group < 75:
        cst_label = 'Orthorhombic'
        cst_marker = 'o'
        cst_marker_size = 3
        cst_zorder = 9
        cst_color = '#0000ff'
        cst_fillstyle = 'none'
    elif space_group < 143:
        cst_label = 'Tetragonal'
        cst_marker = '>'
        cst_marker_size = 4
        cst_zorder = 8
        cst_color = '#00aaaa'
    elif space_group < 168:
        cst_label = 'Trigonal'
        cst_marker = 'v'
        cst_marker_size = 6
        cst_zorder = 8
        cst_color = '#00ff00'
        cst_fillstyle = 'none'
    elif space_group < 195:
        cst_marker = '^'
        cst_label = 'Hexagonal'
        cst_marker_size = 8
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


def main(argv):
    try:
        opts, args = getopt.getopt(argv[1:], "hb:t:o:", ["help", "ylim_bottom=", "ylim_top=", "output_file="])
    except getopt.GetoptError:
        usage(argv[0])
        sys.exit(2)

    # if len(opts) == 0:
    #    usage(argv[0])
    #    sys.exit(2)

    # Default Values
    ylim_bottom = None
    ylim_top = None
    output_file = 'ConvexHull.pdf'

    for opt, arg in opts:
        if opt in ("-h", "--help"):
            usage(argv[0])
            sys.exit()
        elif opt in ("-b", "--ylim_bottom"):
            ylim_bottom = float(arg)
        elif opt in ("-t", "--ylim_top"):
            ylim_top = float(arg)
        elif opt in ("-o", "--output_file"):
            output_file = arg

    if not os.path.isfile('results.json'):
        print('File not found %s' % 'results.json')

    data = json.load(open('results.json'))
    plt.figure(figsize=(11, 8.5))

    energy_left = None
    energy_right = None

    for idata in data:
        if idata['ratio'] == 0:
            if energy_left is None or energy_left > idata['energy_pa']:
                energy_left = idata['energy_pa']
                print('Energy per atom for %2s: %9.3f' % (idata['formula'], energy_left))

        elif idata['ratio'] == 1:
            if energy_right is None or energy_right > idata['energy_pa']:
                energy_right = idata['energy_pa']
                print('Energy per atom for %2s: %9.3f' % (idata['formula'], energy_right))

    if energy_left is None or energy_right is None:
        print('Pure elements not found, formation energy cannot be computed')
        sys.exit(1)

    points = []
    for idata in data:
        x = idata['ratio']
        spcgrp = idata['spcgrp']
        marker, color, lab, m, z, fs = spcgrp_props(spcgrp)
        y = idata['energy_pa'] - (1 - x) * energy_left - x * energy_right
        plt.plot(x, y, marker=marker, ms=m, color=color, fillstyle=fs, zorder=z)
        points.append([x, idata['energy_pa'] - (1 - x) * energy_left - x * energy_right])

    points = np.array(points)
    hull = ConvexHull(points)

    for simplex in hull.simplices:
        if points[simplex, 1][0] <= 0.0 and points[simplex, 1][1] <= 0.0:
            if points[simplex, 1][0] == 0.0 and points[simplex, 1][1] == 0.0:
                continue
            plt.plot(points[simplex, 0], points[simplex, 1], 'k-', zorder=1)

    ylims = plt.ylim()
    if ylim_bottom is None:
        ylim_bottom = ylims[0]
    if ylim_top is None:
        ylim_top = ylims[1]

    for spcgrp in [15, 74, 142, 167, 194, 230]:
        marker, color, lab, m, z, fs = spcgrp_props(spcgrp)
        print('Marker for %12s: %s' % (lab, marker))
        plt.plot(-100, -100, marker, ms=m, fillstyle=fs, color=color, label=lab)

    plt.xlim(-0.05, 1.05)
    print('Limits', ylim_bottom, ylim_top)
    plt.ylim(ylim_bottom, ylim_top)
    plt.legend(loc=9, prop={'size': 10}, numpoints=1)
    plt.subplots_adjust(left=0.12, bottom=0.13, right=0.98, top=0.96, wspace=None, hspace=None)
    plt.xlabel(r'Composition balance')
    plt.ylabel(r'Formation Energy [eV]')
    plt.savefig(output_file)


if __name__ == "__main__":
    main(sys.argv)
