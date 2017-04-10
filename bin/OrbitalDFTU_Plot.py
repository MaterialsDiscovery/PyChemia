#!/usr/bin/env python

from builtins import input
import argparse
import numpy as np
import matplotlib.pyplot as plt
plt.rcdefaults()
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection
import pychemia


def label(xy, text, shiftx):
    px = xy[0] - shiftx
    py = xy[1] - 0.0
    plt.text(px, py, text, va='center', ha="center", family='sans-serif', size=8)

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Orbital DFTU Plotting tool')

    parser.add_argument('-host', type=str, help='Hostname or IP of the Mongo Server', required=True)
    parser.add_argument('-name', type=str, help='Name of Database', required=True)
    parser.add_argument('-dbuser', type=str, help='Username on Database', required=True)
    parser.add_argument('-figname', type=str, help='Path to output figure', default='figure.pdf')
    parser.add_argument('-basepath', type=str, help='Path where calculations are performed', default='.')

    args = parser.parse_args()
    basepath = args.basepath

    passwd = input('Password:')

    db_settings = {'name': args.name, 'host': args.host, 'user': args.dbuser, 'passwd': passwd}
    pcdb = pychemia.db.get_database(db_settings)
    popu = pychemia.population.OrbitalDFTU(pcdb, basepath+'/abinit.in')
    popu.recover()

    ff = pychemia.searcher.FireFly(popu)
    ff.recover()
    print(ff)

    N = ff.current_generation
    M = ff.generation_size

    fig, ax = plt.subplots()
    grid = np.mgrid[0.2:0.8:complex(0, N), 0.2:0.8:complex(0, M)].T

    patches = []
    colors = []

    line_color = {'replace_by_other': 'green', 'replace_by_random': 'cyan', 'promoted': 'black', 'duplicate': 'orange'}
    dx = 0.03

    for j in range(N):
        dup_dx = 0.5*dx
        for i in range(M):
            entry_id = ff.lineage[str(i)][j]
            if ff.population.is_evaluated(entry_id):
                colors.append(ff.population.value(entry_id))
                circle = mpatches.Circle(grid[i, j], 0.3/float(M), ec="none")
                patches.append(circle)
                label(grid[i, j], "%7.2f" % ff.population.value(entry_id), 0.0)

            for ichange in ff.population.pcdb.db.generation_changes.find({'from': entry_id, 'generation': j}):
                if ichange['change'] == 'duplicate':
                    orig = ichange['from']
                    dest = ichange['to']
                    newi = int(ff.lineage_inv[dest])
                    dup_dx += 0.002
                    x, y = np.array([[grid[i, j][0]-1.5*dup_dx, grid[i, j][0]-2*dup_dx,
                                      grid[newi, j][0]-2*dup_dx, grid[newi, j][0]-dx],
                                     [grid[i, j][1], grid[i, j][1], grid[newi, j][1], grid[newi, j][1]]])
                    line = mlines.Line2D(x, y, lw=2., alpha=0.8, color=line_color[ichange['change']],
                                         marker='>', markersize=15, markeredgecolor='none')
                    line.set_markevery([3])
                    ax.add_line(line)
                elif j < N-1:
                    x, y = np.array([[grid[i, j][0]+dx, grid[i, j+1][0]-2*dx], [grid[i, j][1], grid[i, j+1][1]]])
                    line = mlines.Line2D(x, y, lw=5., alpha=0.3, color=line_color[ichange['change']])
                    label(0.5*(grid[i, j]+grid[i, j+1]), ichange['change'], 0.0)
                    ax.add_line(line)

    collection = PatchCollection(patches, cmap=plt.cm.hsv, alpha=0.3)
    collection.set_array(np.array(colors))
    ax.add_collection(collection)

    plt.subplots_adjust(left=0, right=1, bottom=0, top=1)
    plt.axis('equal')
    plt.axis('off')
#    plt.show()
    plt.savefig(args.figname)
