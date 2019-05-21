#!/usr/bin/env python

import os
import argparse
import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
plt.rcdefaults()
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection
import pychemia
import bson


def label(xy, text, shiftx):
    px = xy[0] - shiftx
    py = xy[1] - 0.0
    plt.text(px, py, text, va='center', ha="center", family='sans-serif', size=8)


def compare_params(path):

    abinitout = popu.get_final_abinit_out(path)

    if not os.path.isfile(path+os.sep+'abinit.in'):
        print('ERROR: No abinit.in found at %s' % path)
        return

    # For making easier to see the values
    np.set_printoptions(linewidth=200, suppress=True)

    # Reading the INPUT
    abi = pychemia.code.abinit.AbinitInput(path+os.sep+'abinit.in')
    idmatpawu = np.array(abi['dmatpawu']).reshape(-1, 5, 5)
    iparams = pychemia.population.orbitaldftu.dmatpawu2params(idmatpawu, 5)

    # Reading the OUTPUT
    dmatpawu = pychemia.population.orbitaldftu.get_final_dmatpawu(abinitout)
    odmatpawu = np.array(dmatpawu).reshape(-1, 5, 5)
    oparams = pychemia.population.orbitaldftu.dmatpawu2params(odmatpawu, 5)

    print('PARAMETRIC REPRESENTATION found at %s' % abinitout)
    for i in sorted(list(iparams.keys())):
        print(i)
        print('input')
        print(iparams[i])
        print('output')
        print(i)
        print(oparams[i])

    abo = pychemia.code.abinit.AbinitOutput(abinitout)
    if not abo.is_finished:
        print('This output is not finished')
    try:
        nres2 = abo.get_energetics()['nres2'][-1]
        etot = abo.get_energetics()['etot'][-1]
        nscf = len(abo.get_energetics()['etot'])
        print("%30s ETOT: %15.6f NRES2: %15.6e NumSCF: %3d" % (path, etot, nres2, nscf))
    except:
        print("ERROR: Could not get energetics from %s" % abinitout)


def create_population():
    popu.random_population(64)


def prepare_folders(scrpath):
    for i in popu.members:
        popu.prepare_folder(i, workdir=scrpath, source_dir=scrpath)


def check_status(basepath, dirs):

    print("%-40s %15s %15s %4s" % ("ABINIT output", "ETOT", 'nres2', 'nSCF'))
    for i in dirs:
        path = basepath+os.sep+i
        abinitout = popu.get_final_abinit_out(path)
        if abinitout is None:
            continue

        abo = pychemia.code.abinit.AbinitOutput(abinitout)
        if not abo.is_finished:
            continue
        try:
            nres2 = abo.get_energetics()['nres2'][-1]
            etot = abo.get_energetics()['etot'][-1]
            nscf = len(abo.get_energetics()['etot'])
            print("%-40s %15.6f %15.6e %4d" % (abinitout, etot, nres2, nscf))
        except:
            print("ERROR: Could not get final energetics from %s" % abinitout)


def plot_status(basepath, dirs, nstep=50):
    print('Plotting Status...')
    X = {}
    Y = {}

    # Get the data:
    for path in dirs:
        X[path] = np.array([])
        Y[path] = {}
        Y[path]['nres2'] = np.array([])
        Y[path]['etot'] = np.array([])
        Y[path]['delta'] = np.array([])

        for i in range(10):
            if os.path.isfile(basepath+os.sep+path+os.sep+'abinit_%d.out' % i):
                abo = pychemia.code.abinit.AbinitOutput(basepath+os.sep+path+os.sep+'abinit_%d.out' % i)
                if not abo.is_finished:
                    print('This output is not finished')
                try:
                    energ = abo.get_energetics()
                    nres2 = energ['nres2'][-1]
                    etot = energ['etot'][-1]
                    nscf = len(energ['etot'])
                    x = np.arange(nstep*i, nstep*i + len(energ['nres2']))
                    yetot = np.array(energ['etot'])
                    ynres2 = np.array(energ['nres2'])
                    ydelta = np.array(np.abs(energ['deltaEh']))
                    X[path] = np.concatenate((X[path], x))
                    Y[path]['etot'] = np.concatenate((Y[path]['etot'], yetot))
                    Y[path]['nres2'] = np.concatenate((Y[path]['nres2'], ynres2))
                    Y[path]['delta'] = np.concatenate((Y[path]['delta'], ydelta))
                    print("%s ETOT:%15.6f NRES2=%15.6e Num SCF=%3d" % (path, etot, nres2, nscf))
                except:
                    print("%s Failed processing output" % path)

    # RESIDUAL
    plt.figure(figsize=(8, 11))
    plt.subplots_adjust(left=0.1, bottom=0.05, right=0.99, top=0.99,
                        wspace=None, hspace=None)
    for path in dirs:
        if len(X) > 0:
            plt.semilogy(X[path], Y[path]['nres2'], label=path[-4:])

    for i in nstep*np.arange(10):
        plt.semilogy([i, i], [1E-19, 1E5], '0.5')

    plt.ylim(1E-19, 1E5)
    plt.xlim(0, max([max(X[path]) for path in dirs]))
    plt.legend()
    plt.xlabel("SCF iteration")
    plt.ylabel("Density Residual$^2$")
    plt.savefig('Residual.pdf')

    # ETOT
    plt.figure(figsize=(8, 11))
    for path in dirs:
        if len(X) > 0:
            plt.plot(X[path], Y[path]['etot'], label=path[-4:])

    plt.subplots_adjust(left=0.1, bottom=0.05, right=0.99, top=0.99,
                        wspace=None, hspace=None)

    for i in nstep*np.arange(10):
        plt.plot([i, i], [-992, -990], '0.5')

    plt.ylim(-992, -990)
    plt.xlim(0, max([max(X[path]) for path in dirs]))
    plt.legend()
    plt.xlabel("SCF iteration")
    plt.ylabel("Energy")
    plt.savefig('ETotal.pdf')

    # deltaEh
    plt.figure(figsize=(8, 11))
    for path in dirs:
        if len(X) > 0:
            plt.semilogy(X[path], Y[path]['delta'], label=path[-4:])

    plt.subplots_adjust(left=0.1, bottom=0.05, right=0.99, top=0.99,
                        wspace=None, hspace=None)

    for i in nstep*np.arange(10):
        plt.semilogy([i, i], [1E-13, 1E3], '0.5')

    plt.ylim(1E-13, 1E3)
    plt.xlim(0, max([max(X[path]) for path in dirs]))
    plt.legend()
    plt.xlabel("SCF iteration")
    plt.ylabel("Delta Energy")
    plt.savefig('deltaEh.pdf')


def plot_polar():
    print('Computing duplicates in %d candidates' % len(popu.evaluated))
    dupes = popu.get_duplicates(popu.evaluated)
    print('Number of duplicates    : %d' % len(dupes))
    nodupes = [x for x in popu.evaluated if x not in dupes]
    print('Number of not duplicates: %d' % len(nodupes))

    if candidates is not None:
        nodupes = candidates
    
    fig = plt.figure(figsize=(20, 2.01*len(nodupes)))

    etots = np.array([popu.value(x) for x in nodupes])
    sort_entries = np.array(nodupes)[etots.argsort()]

    N = len(nodupes)
    M = 10

    grid = np.mgrid[0.02:0.87:complex(0, M), 0.02:0.84:complex(0, N)].T
    print("List of duplicates to plot:")
    for i in nodupes:
        print(i)

    for i in range(N):

        entry_id = sort_entries[i]
        pm = popu.get_correlation_params(entry_id)
        etot = popu.value(entry_id)
        ea = np.array(pm['euler_angles'])
        nres2 = popu.get_entry(entry_id, {'properties.nres2': 1})['properties']['nres2']

        fig.text(0.03, grid[i, 0][1]+0.08, "%7.2f" % etot, va='center', ha='center', fontsize='24', rotation='vertical',
                 color='red')
        fig.text(0.01, grid[i, 0][1]+0.08, "%9.2E" % nres2, va='center', ha='center', fontsize='16', rotation='vertical',
                 color='blue')

        radii = np.array([1, 1, 1, 1])
        #colors = np.array([0, 1, 2, 3])
        colors = ['red', 'blue', 'fuchsia', 'yellow']
        width = np.array([0.15, 0.15, 0.1, 0.1])

        for j in range(M):

            theta = ea[:, j]

            ax = fig.add_axes([grid[i, j][0], grid[i, j][1], 0.146, 0.146], projection='polar')
            ax.set_axis_bgcolor('mintcream')
            #ax.yaxis.set_tick_params(labelsize=0)
            #ax.xaxis.set_tick_params(labelsize=0)
            ax.get_xaxis().set_visible(False)
            ax.get_yaxis().set_visible(False)
            ax.spines['polar'].set_visible(True)
            bars = ax.bar(theta, radii, width=width, bottom=0.0, linewidth=0)

            # Use custom colors and opacity
            for r, bar in zip(colors, bars):
                #bar.set_facecolor(plt.cm.Paired(r))
                bar.set_facecolor(r)
                bar.set_alpha(0.9)

    plt.subplots_adjust(left=0, right=1, bottom=0, top=1)
    #plt.axis('equal')
    #plt.axis('off')
    plt.savefig(figname+'_polar.png')


def plot_evolution():
    N = ff.current_generation
    M = ff.generation_size

    fig, ax = plt.subplots(figsize=(20, 20))
    grid = np.mgrid[0.1:0.9:complex(0, N), 0.1:0.9:complex(0, M)].T

    patches = []
    colors = []

    line_color = {'replace_by_other': 'green', 'replace_by_random': 'cyan', 'promoted': 'black', 'duplicate': 'orange'}
    dx = 0.01

    for j in range(N):
        dup_dx = 0.5 * dx
        for i in range(M):
            entry_id = ff.lineage[str(i)][j]

            nres2 = ff.population.get_entry(entry_id, {'properties.nres2': 1})['properties']['nres2']
            if nres2 < 1E14:
                lw = 1
                ls = 'solid'
            else:
                lw = 1
                ls = 'dotted'

            if ff.population.is_evaluated(entry_id):
                colors.append(ff.population.value(entry_id))
                circle = mpatches.Circle(grid[i, j], 0.4 / float(M), ec="black", linewidth=lw, linestyle=ls)
                patches.append(circle)
                label(grid[i, j], "%7.2f" % ff.population.value(entry_id), 0.0)

            for ichange in ff.population.pcdb.db.generation_changes.find({'from': entry_id, 'generation': j}):
                if ichange['change'] == 'duplicate':
                    orig = ichange['from']
                    dest = ichange['to']
                    newi = int(ff.lineage_inv[dest])
                    dup_dx += dx/10.0
                    x, y = np.array([[grid[i, j][0] - 1.5 * dup_dx, grid[i, j][0] - 2 * dup_dx,
                                      grid[newi, j][0] - 2 * dup_dx, grid[newi, j][0] - dx],
                                     [grid[i, j][1], grid[i, j][1], grid[newi, j][1], grid[newi, j][1]]])
                    line = mlines.Line2D(x, y, lw=1., alpha=0.8, color=line_color[ichange['change']],
                                         marker='>', markersize=5, markeredgecolor='none')
                    line.set_markevery([3])
                    ax.add_line(line)
                elif j < N - 1:
                    x, y = np.array(
                        [[grid[i, j][0] + dx, grid[i, j + 1][0] - 2 * dx], [grid[i, j][1], grid[i, j + 1][1]]])
                    line = mlines.Line2D(x, y, lw=5., alpha=0.3, color=line_color[ichange['change']])
                    # label(0.5*(grid[i, j]+grid[i, j+1]), ichange['change'], 0.0)
                    ax.add_line(line)

    collection = PatchCollection(patches, cmap=plt.cm.hsv, alpha=0.3)
    collection.set_array(np.array(colors))
    ax.add_collection(collection)

    plt.subplots_adjust(left=0, right=1, bottom=0, top=1)
    plt.axis('equal')
    plt.axis('off')
    plt.savefig(figname+'_evo.pdf')


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Orbital DFTU Plotting tool')

    parser.add_argument('-host', type=str, help='Hostname or IP of the Mongo Server', required=True)
    parser.add_argument('-name', type=str, help='Name of Database', required=True)
    parser.add_argument('-dbuser', type=str, help='Username on Database', required=True)
    parser.add_argument('-figname', type=str, help='Path to output figure', default=None)
    parser.add_argument('-basepath', type=str, help='Path where calculations are performed', default='.')
    parser.add_argument('-candidates', type=str, help='IDs for candidates', nargs='+')

    args = parser.parse_args()
    basepath = args.basepath

    if args.figname is None:
        figname = args.name
    else:
        figname = args.name

    if figname[-4] == '.':
        figname = figname[:-4]

    passwd = input('Password:')

    db_settings = {'name': args.name, 'host': args.host, 'user': args.dbuser, 'passwd': passwd}
    pcdb = pychemia.db.get_database(db_settings)
    popu = pychemia.population.OrbitalDFTU(pcdb, basepath+'/abinit.in')
    popu.recover()

    ff = pychemia.searcher.FireFly(popu)
    ff.recover()
    print(ff)

    if len(args.candidates) > 0:
        candidates = [bson.ObjectId(x) for x in args.candidates]
    else:
        candidates = None

    for i in candidates:
        assert(popu.get_entry(i) is not None)
        entry = popu.get_entry(i)
        print("%s %f" % (entry['_id'], entry['properties']['etot']))
        
    print('Plotting Evolution')
    plot_evolution()
    print('Plotting Euler Angles')
    plot_polar()
