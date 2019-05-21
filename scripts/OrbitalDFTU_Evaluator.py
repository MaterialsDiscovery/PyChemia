#!/usr/bin/env python

import os
import sys
import json
import argparse
import numpy as np
import matplotlib
import multiprocessing
import logging

matplotlib.use('agg')
import matplotlib.pyplot as plt
from pychemia.code.abinit import AbinitInput, AbinitOutput
from pychemia.population.orbitaldftu import dmatpawu2params, params2dmatpawu, OrbitalDFTU
from pychemia.population.orbitaldftu import get_final_abinit_out
from pychemia.searcher import FireFly
from pychemia.db import get_database


def compare_params(path):
    if not os.path.isfile(path + os.sep + 'abinit.in'):
        print('ERROR: No abinit.in found at %s' % path)
        return
    # For making easier to see the values
    np.set_printoptions(linewidth=200, suppress=True)

    # Reading the INPUT
    abi = AbinitInput(path + os.sep + 'abinit.in')

    if 'lpawu' not in abi.variables:
        raise ValueError("Variable lpawu not found")
    ndim = 2 * max(abi['lpawu']) + 1

    idmatpawu = np.array(abi['dmatpawu']).reshape(-1, ndim, ndim)
    iparams = dmatpawu2params(idmatpawu, ndim)

    # Reading the OUTPUT
    abinitout = get_final_abinit_out(path)
    abo = AbinitOutput(abinitout)
    dmatpawu = abo.get_final_dmatpawu()
    odmatpawu = np.array(dmatpawu).reshape(-1, ndim, ndim)
    oparams = dmatpawu2params(odmatpawu, ndim)

    print('DMATPAWU')
    print('input')
    print(idmatpawu)
    print('output')
    print(odmatpawu)

    print('PARAMETRIC REPRESENTATION found at %s' % abinitout)
    for i in sorted(list(iparams.keys())):
        print(i)
        print('input')
        print(iparams[i])
        print('output')
        print(i)
        print(oparams[i])

    abo = AbinitOutput(abinitout)
    if not abo.is_finished:
        print('This output is not finished')
    try:
        nres2 = abo.get_energetics()['nres2'][-1]
        etot = abo.get_energetics()['etot'][-1]
        nscf = len(abo.get_energetics()['etot'])
        print("%30s ETOT: %15.6f NRES2: %15.6e NumSCF: %3d" % (path, etot, nres2, nscf))
    except:
        print("ERROR: Could not get energetics from %s" % abinitout)


def check_status(basepath):
    dirs = [x for x in os.listdir(basepath)
            if os.path.isdir(basepath + os.sep + x) and os.path.isfile(basepath + os.sep + x + os.sep + 'abinit.in')]
    print("%-40s %15s %15s %4s" % ("ABINIT output", "ETOT", 'nres2', 'nSCF'))
    for i in dirs:
        path = basepath + os.sep + i
        abinitout = get_final_abinit_out(path)
        if abinitout is None:
            continue

        abo = AbinitOutput(abinitout)
        if not abo.is_finished:
            continue
        try:
            nres2 = abo.get_energetics()['nres2'][-1]
            etot = abo.get_energetics()['etot'][-1]
            nscf = len(abo.get_energetics()['etot'])
            print("%-40s %15.6f %15.6e %4d" % (abinitout, etot, nres2, nscf))
        except:
            print("ERROR: Could not get final energetics from %s" % abinitout)


def plot_polar(popu, basepath):
    print('Plotting Euler Angles...')

    dirs = [x for x in os.listdir(basepath)
            if os.path.isdir(basepath + os.sep + x) and os.path.isfile(basepath + os.sep + x + os.sep + 'abinit.in')]
    fig = plt.figure(figsize=(21, 1.2 * len(dirs)))
    plt.subplots_adjust(left=0.0, bottom=0.0, right=0.99, top=0.99, wspace=None, hspace=None)
    etots = []

    for idir in dirs:
        pm, etot = popu.get_final_properties(basepath + os.sep + idir)
        etots.append(etot)

    etots = np.array(etots)
    sort_dirs = np.array(dirs)[etots.argsort()]

    index = 0.0
    for idir in sort_dirs:

        pm, etot = popu.get_final_properties(basepath + os.sep + idir)
        ea = np.array(pm['final_dmat']['euler_angles'])
        # print(idir,etot)
        # print(ea.shape)
        nangles = ea.shape[1]

        for j in range(nangles):
            theta = ea[:, j]
            dim = len(theta)
            radii = np.ones(dim)
            colors = np.arange(dim)
            width = 0.1 * np.ones(dim)

            ax = fig.add_axes([float(j) / nangles, index / (len(dirs) + 1), 1.0 / nangles, 1.0 / nangles],
                              projection='polar')
            ax.yaxis.set_tick_params(labelsize=0)
            ax.xaxis.set_tick_params(labelsize=0)
            ax.get_xaxis().set_visible(False)
            ax.get_yaxis().set_visible(False)
            ax.spines['polar'].set_visible(True)
            bars = ax.bar(theta, radii, width=width, bottom=0.0)

            # Use custom colors and opacity
            for r, bar in zip(colors, bars):
                bar.set_facecolor(plt.cm.viridis(float(r) / nangles))
                bar.set_alpha(0.9)

        index += 1.0

    plt.savefig('OrbitalPolar.pdf')
    plt.show()


def plot_status(basepath):
    dirs = [x for x in os.listdir(basepath)
            if os.path.isdir(basepath + os.sep + x) and os.path.isfile(basepath + os.sep + x + os.sep + 'abinit.in')]

    abi = AbinitInput(basepath + os.sep + 'abinit.in')
    nstep = abi['nstep']

    print('Plotting Status...')
    xx = {}
    yy = {}

    # Get the data:
    max_nruns = 0
    for path in dirs:
        xx[path] = np.array([])
        yy[path] = {}
        yy[path]['nres2'] = np.array([])
        yy[path]['etot'] = np.array([])
        yy[path]['delta'] = np.array([])

        for i in range(100):
            if os.path.isfile(basepath + os.sep + path + os.sep + 'abinit_%02d.txt' % i):
                abo = AbinitOutput(basepath + os.sep + path + os.sep + 'abinit_%02d.txt' % i)
                if not abo.is_finished:
                    print('This output is not finished')
                    continue
                if max_nruns < i:
                    max_nruns = i
                try:
                    energ = abo.get_energetics()
                except:
                    raise RuntimeError(
                        "Failed procesing: %s" % (basepath + os.sep + path + os.sep + 'abinit_%02d.txt' % i))

                nres2 = energ['nres2'][-1]
                etot = energ['etot'][-1]
                nscf = len(energ['etot'])
                x = np.arange(nstep * i, nstep * i + len(energ['nres2']))
                yetot = np.array(energ['etot'])
                ynres2 = np.array(energ['nres2'])
                ydelta = np.array(np.abs(energ['deltaEh']))
                xx[path] = np.concatenate((xx[path], x))
                yy[path]['etot'] = np.concatenate((yy[path]['etot'], yetot))
                yy[path]['nres2'] = np.concatenate((yy[path]['nres2'], ynres2))
                yy[path]['delta'] = np.concatenate((yy[path]['delta'], ydelta))
                print("%s ETOT:%15.6f NRES2=%15.6e Num SCF=%3d" % (path, etot, nres2, nscf))

    # RESIDUAL
    plt.figure(figsize=(8, 11))
    miny = 1E-1
    maxy = 1E-16
    plt.subplots_adjust(left=0.1, bottom=0.05, right=0.99, top=0.99,
                        wspace=None, hspace=None)
    for path in dirs:
        if len(xx) > 0:
            plt.semilogy(xx[path], yy[path]['nres2'], label=path[-4:], lw=0.1)
            if miny > min(yy[path]['nres2']):
                miny = min(yy[path]['nres2'])
            if maxy < max(yy[path]['nres2']):
                maxy = max(yy[path]['nres2'])

    plt.ylim(miny, maxy)

    for i in nstep * np.arange(max_nruns + 1):
        plt.semilogy([i, i], [miny, maxy], '0.5', lw=0.1)

    plt.xlim(0, max([max(xx[path]) for path in dirs]))
    plt.legend()
    plt.xlabel("SCF iteration")
    plt.ylabel("Density Residual$^2$")
    plt.savefig('Orbital_Residual.pdf')

    # ETOT
    miny = 1E6
    maxy = -1E6
    avg = 0
    minlen = 100
    plt.figure(figsize=(8, 11))
    for path in dirs:
        if len(xx) > 0:
            plt.plot(xx[path], yy[path]['etot'], label=path[-4:])
            if len(yy[path]['etot']) < minlen:
                minlen = len(yy[path]['etot'])
                avg = np.average(yy[path]['etot'][-int(minlen / 2):])
            if miny > min(yy[path]['etot'][-int(minlen / 2):]):
                miny = min(yy[path]['etot'][-int(minlen / 2):])
            if maxy < max(yy[path]['etot'][-int(minlen / 2):]):
                maxy = max(yy[path]['etot'][-int(minlen / 2):])

    plt.subplots_adjust(left=0.15, bottom=0.05, right=0.95, top=0.95,
                        wspace=None, hspace=None)

    newminy = miny - 0.1 * (maxy - miny)
    newmaxy = maxy + 0.1 * (maxy - miny)
    miny = newminy
    maxy = newmaxy

    plt.ylim(miny, maxy)
    for i in nstep * np.arange(max_nruns + 1):
        plt.plot([i, i], [miny, maxy], '0.5', lw=0.1)

    plt.xlim(0, max([max(xx[path]) for path in dirs]))
    plt.legend()
    plt.xlabel("SCF iteration")
    plt.ylabel("Energy")
    plt.savefig('Orbital_ETotal.pdf')

    # deltaEh
    plt.figure(figsize=(8, 11))
    miny = 1E-1
    for path in dirs:
        if len(xx) > 0:
            plt.semilogy(xx[path], yy[path]['delta'], label=path[-4:], lw=0.1)
            if miny > min(yy[path]['delta']):
                miny = min(yy[path]['delta'])

    plt.subplots_adjust(left=0.1, bottom=0.05, right=0.99, top=0.99,
                        wspace=None, hspace=None)

    for i in nstep * np.arange(max_nruns + 1):
        plt.semilogy([i, i], [miny, 1E3], '0.5', lw=0.1)

    plt.ylim(miny, 1E3)
    plt.xlim(0, max([max(xx[path]) for path in dirs]))
    plt.legend()
    plt.xlabel("SCF iteration")
    plt.ylabel("Delta Energy")
    plt.savefig('Orbital_deltaEh.pdf')


def create_population(num_candidates):
    popu.random_population(num_candidates)


def safe_read_json(filename):
    """
    Safely read a given filename, extract and returns the associated dictionary
    from the JSON file
    """
    if not os.path.exists(filename):
        raise ValueError("ERROR: Could not read file: %s" % filename)
    rf = open(filename)
    try:
        data = json.load(rf)
    except ValueError:
        raise ValueError("ERROR: File is not in proper JSON format: %s" % filename)
    rf.close()
    return data


def prepare_folders(scrpath):
    for i in popu.members:
        popu.prepare_folder(i, workdir=scrpath, source_dir=scrpath)


def set_populations(parser):

    return parser.add_argument('-p', type=str, nargs='+', required=True, metavar='JSON_file',
                               help='Population settings (JSON file), includes settings to connect to server and the '
                                    'population')


def evaluate(db_settings, queue_settings, abipath):
    pcdb = get_database(db_settings)
    print("[%s] Path for abinit.in: %s" % (pcdb.name, abipath))
    popu = OrbitalDFTU(pcdb, abipath + os.sep + 'abinit.in')
    popu.evaluator(queue_settings, abipath)


def search(db_settings, search_settings, abipath):
    pcdb = get_database(db_settings)
    print("[%s] Path for abinit.in: %s" % (pcdb.name, abipath))
    popu = OrbitalDFTU(pcdb, abipath + os.sep + 'abinit.in')

    if 'generation_size' in search_settings:
        generation_size = search_settings.pop('generation_size')
    else:
        generation_size = 32

    if 'stabilization_limit' in search_settings:
        stabilization_limit = search_settings.pop('stabilization_limit')
    else:
        stabilization_limit = 10

    fire = FireFly(popu, params=search_settings, generation_size=generation_size,
                   stabilization_limit=stabilization_limit)
    fire.run()


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Orbital DFTU Evaluator and Analysis Tool')

    subparsers = parser.add_subparsers(help='commands', dest='subparser_name')

    # The create command
    create_parser = subparsers.add_parser('create', help='Create the database')

    set_populations(create_parser)

    # The populate command
    populate_parser = subparsers.add_parser('populate', help='Add candidates to the population (used for testing)')

    set_populations(populate_parser)

    populate_parser.add_argument('-clean', action='store_true',
                                 help='If true, clean the entire database before populate',
                                 required=False, default=False)

    populate_parser.add_argument('-size', type=int,
                                 help='Number of candidates to populate (default: 16)',
                                 required=False, default=16)
    populate_parser.add_argument('-basepath', type=str,
                                 help='Path where calculations are performed',
                                 required=False, default='.')

    # A run command
    run_parser = subparsers.add_parser('run', help='Run Evaluator')

    set_populations(run_parser)

    run_parser.add_argument('-queue_settings', type=str,
                            help='Filename with PBS settings for launching jobs',
                            required=False, default='queue.json')
    run_parser.add_argument('-basepath', type=str,
                            help='Path where calculations are performed (default: current folder)',
                            required=False, default='.')

    # A searcher command
    searcher_parser = subparsers.add_parser('search', help='Run PyChemia Global Searcher')

    set_populations(searcher_parser)

    searcher_parser.add_argument('-search_settings', type=str,
                                 help='Filename with PBS settings for launching jobs',
                                 required=False, default='searcher.json')
    searcher_parser.add_argument('-basepath', type=str,
                                 help='Path where calculations are performed (default: current folder)',
                                 required=False, default='.')

    # The plot command
    plot_parser = subparsers.add_parser('plot', help='Generate several plots')

    set_populations(plot_parser)

    plot_parser.add_argument('-basepath', type=str,
                             help='Path where calculations are performed',
                             required=False, default='.')

    args = parser.parse_args()
    print(args)

    # check all settings in args.p
    dbs = []
    for dbi_file in args.p:
        dbi = safe_read_json(dbi_file)
        print("DB settings: %s" % dbi)
        assert('name' in dbi)
        assert('u' in dbi)
        assert('j' in dbi)
        dbs.append(dbi)

    if args.subparser_name == 'create':
        
        pcdbs = []
        for dbi in dbs:
            pcdb = get_database(dbi)
            pcdbs.append(pcdb)
            print(pcdb)
            print(pcdb.entries.count())

    if args.subparser_name == 'run':
        queue_settings = safe_read_json(args.queue_settings)
        print("PBS settings: %s" % queue_settings)

    # General actions for 'populate', 'run', 'search' and 'plot'
    if args.subparser_name in ['run', 'plot', 'populate', 'search']:

        if not os.path.isdir(args.basepath) or not os.path.isfile(args.basepath + os.sep + 'abinit.in'):
            print('ERROR: Wrong basepath %s, directory must exist and contain a abinit.in file' % args.basepath)
            parser.print_help()
            sys.exit(1)

        popu = {}
        for dbi in dbs:
            name = dbi['name']
            pcdb = get_database(dbi)
            if not os.path.isdir(args.basepath + os.sep + name):
                os.mkdir(args.basepath + os.sep + name)
            abi = AbinitInput(args.basepath + os.sep + 'abinit.in')
            
            abi['upawu'] = ""
            for i in dbi['u']:
                abi['upawu'] += str(i) + " "
            abi['upawu'] += 'eV'
                
            abi['jpawu'] = ""
            for i in dbi['j']:
                abi['jpawu'] += str(i) + " "
            abi['jpawu'] += 'eV'

            abipath = args.basepath + os.sep + name + os.sep + 'abinit.in'
            abi.write(abipath)

            abifiles = args.basepath + os.sep + name + os.sep + 'abinit.files'
            if os.path.lexists(abifiles):
                os.remove(abifiles)
            os.symlink(os.path.abspath(args.basepath + os.sep + 'abinit.files'), abifiles)
            for i in [x for x in os.listdir(args.basepath) if x[-3:] == 'xml']:
                psp = args.basepath + os.sep + name + os.sep + i
                if os.path.lexists(psp):
                    os.remove(psp)
                os.symlink(os.path.abspath(args.basepath + os.sep + i), psp)

            popu[name] = OrbitalDFTU(pcdb, abipath)

    # Specific actions for 'populate', 'run', 'search' and 'plot'
    if args.subparser_name == 'populate':

        for dbi in dbs:
            name = dbi['name']
            if args.clean:
                popu[name].clean()
            popu[name].random_population(args.size)
            print("[%s] Total number of candidates: %d" % (name, len(popu[name])))

    elif args.subparser_name == 'run':
        sprocs = {}
        for dbi in dbs:
            name = dbi['name']
            basepath = args.basepath + os.sep + name

            if os.path.exists(basepath + os.sep + 'ERROR'):
                print('ERROR: Something was wrong with %s' % abipath)
            else:
                sprocs[name] = multiprocessing.Process(target=evaluate, args=(dbi, queue_settings, basepath))
                sprocs[name].start()

        for dbi in dbs:            
            name = dbi['name']
            sprocs[name].join()

    elif args.subparser_name == 'search':
        logging.basicConfig(level=logging.DEBUG)
        search_settings = safe_read_json(args.search_settings)
        print("Search settings from file: %s \n%s" % (args.search_settings, search_settings))

        sprocs = {}
        for dbi in dbs:
            name = dbi['name']
            basepath = args.basepath + os.sep + name
            sprocs[name] = multiprocessing.Process(target=search, args=(dbi, search_settings, basepath))
            sprocs[name].start()

        print(list(sprocs.keys()))

        for dbi in dbs:
            name = dbi['name']
            if name not in sprocs:
                print('ERRROR: %s' % str(sprocs))
            sprocs[name].join()

    elif args.subparser_name == 'plot':
        check_status(args.basepath)
        plot_status(args.basepath)
        plot_polar(popu, args.basepath)
