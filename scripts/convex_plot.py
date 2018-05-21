#!/usr/bin/env python

import os
import sys
import json
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull
import argparse


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


def formation_energy(energy, x, energy_left, energy_right):
    return energy - (1 - x) * energy_left - x * energy_right


def create_convex(bottom, top, energy_left, energy_right, input, output):

    if not os.path.isfile(input):
        print('File not found %s' % input)

    data = json.load(open(input))
    plt.figure(figsize=(11, 8.5))

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
        y = formation_energy(idata['energy_pa'], x, energy_left, energy_right)
        plt.plot(x, y, marker=marker, ms=m, color=color, fillstyle=fs, zorder=z)
        points.append([x, formation_energy(idata['energy_pa'], x, energy_left, energy_right)])

    points.append([0.0, 0.0])
    points.append([1.0, 0.0])

    points = np.array(points)
    hull = ConvexHull(points)

    for simplex in hull.simplices:
        if points[simplex, 1][0] <= 0.0 and points[simplex, 1][1] <= 0.0:
            if points[simplex, 1][0] == 0.0 and points[simplex, 1][1] == 0.0:
                continue
            plt.plot(points[simplex, 0], points[simplex, 1], 'k-', zorder=1)

    ylims = plt.ylim()
    if bottom is None:
        bottom = ylims[0]
    if top is None:
        top = ylims[1]

    for spcgrp in [15, 74, 142, 167, 194, 230]:
        marker, color, lab, m, z, fs = spcgrp_props(spcgrp)
        print('Marker for %12s: %s' % (lab, marker))
        plt.plot(-100, -100, marker, ms=m, fillstyle=fs, color=color, label=lab)

    plt.xlim(-0.05, 1.05)
    print('Limits', bottom, top)
    plt.ylim(bottom, top)
    plt.legend(loc=9, prop={'size': 10}, numpoints=1)
    plt.subplots_adjust(left=0.12, bottom=0.13, right=0.98, top=0.96, wspace=None, hspace=None)
    plt.xlabel(r'Composition balance')
    plt.ylabel(r'Formation Energy [eV]')
    plt.savefig(output)
    return plt.gcf()


if __name__ == "__main__":

    description = """Collect structures from several databases for plotting convex hulls"""

    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-t', '--top',
                        default=None, metavar='ymax', type=float,
                        help='Maximum value of formation energy')
    parser.add_argument('-b', '--bottom',
                        default=None, metavar='ymax', type=float,
                        help='Minimum value of formation energy')
    parser.add_argument('-l', '--left_energy_pa',
                        default=None, metavar='energy_pa', type=float,
                        help='Energy per atom of left specie')
    parser.add_argument('-r', '--right_energy_pa',
                        default=None, metavar='energy_pa', type=float,
                        help='Energy per atom right specie')
    parser.add_argument('-i', '--input',
                        default='convex.json', metavar='convex.json', type=str,
                        help='Input file for the Convex Hull (JSON file)')
    parser.add_argument('-o', '--output',
                        default='convex.pdf', metavar='convex.pdf', type=str,
                        help='Output file (default: convex.pdf)')

    args = parser.parse_args()
    if not os.path.isfile(args.input):
        parser.print_help()
        exit(1)
    print(args)

    create_convex(args.bottom, args.top, args.left_energy_pa, args.right_energy_pa, args.input, args.output)
