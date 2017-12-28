#!/usr/bin/env python

from builtins import input
import os
import sys
import json
import argparse
import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from pychemia.code.abinit import AbinitInput, AbinitOutput
from pychemia.population.orbitaldftu import dmatpawu2params, params2dmatpawu, OrbitalDFTU
from pychemia.population.orbitaldftu import get_final_abinit_out
from pychemia.db import get_database

def compare_params(path):

    if not os.path.isfile(path+os.sep+'abinit.in'):
        print('ERROR: No abinit.in found at %s' % path)
        return
    # For making easier to see the values
    np.set_printoptions(linewidth=200, suppress=True)

    # Reading the INPUT
    abi=AbinitInput(path+os.sep+'abinit.in')

    if 'lpawu' not in abi.variables:
        raise ValueError("Variable lpawu not found")
    ndim = 2*max(abi['lpawu'])+1

    idmatpawu=np.array(abi['dmatpawu']).reshape(-1,ndim,ndim)
    iparams=dmatpawu2params(idmatpawu,ndim)

    # Reading the OUTPUT
    abinitout = get_final_abinit_out(path)
    abo=AbinitOutput(abinitout)
    dmatpawu= abo.get_final_dmatpawu()
    odmatpawu=np.array(dmatpawu).reshape(-1,ndim,ndim)
    oparams=dmatpawu2params(odmatpawu,ndim)

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

    abo=AbinitOutput(abinitout)
    if not abo.is_finished:
        print('This output is not finished')
    try:
        nres2=abo.get_energetics()['nres2'][-1]
        etot=abo.get_energetics()['etot'][-1]
        nscf=len(abo.get_energetics()['etot'])
        print("%30s ETOT: %15.6f NRES2: %15.6e NumSCF: %3d" % (path,etot, nres2, nscf))
    except:
        print("ERROR: Could not get energetics from %s" % abinitout)

def check_status(basepath):

    dirs=[ x for x in os.listdir(basepath) 
           if os.path.isdir(basepath+os.sep+x) and os.path.isfile(basepath+os.sep+x+os.sep+'abinit.in') ]
    print("%-40s %15s %15s %4s" % ("ABINIT output", "ETOT", 'nres2', 'nSCF'))
    for i in dirs:
        path=basepath+os.sep+i
        abinitout = get_final_abinit_out(path)
        if abinitout is None:
            continue

        abo=AbinitOutput(abinitout)
        if not abo.is_finished:
            continue
        try:
            nres2=abo.get_energetics()['nres2'][-1]
            etot=abo.get_energetics()['etot'][-1]
            nscf=len(abo.get_energetics()['etot'])
            print("%-40s %15.6f %15.6e %4d" % (abinitout, etot, nres2, nscf))
        except:
            print("ERROR: Could not get final energetics from %s" % (abinitout))


def plot_polar(popu, basepath):
    print('Plotting Euler Angles...')

    dirs=[ x for x in os.listdir(basepath) 
           if os.path.isdir(basepath+os.sep+x) and os.path.isfile(basepath+os.sep+x+os.sep+'abinit.in') ]
    fig=plt.figure(figsize=(21,1.2*len(dirs)))
    plt.subplots_adjust(left=0.0, bottom=0.0, right=0.99, top=0.99, wspace=None, hspace=None)
    etots=[]

    for idir in dirs:
        pm, etot=popu.get_final_properties(basepath+os.sep+idir)
        etots.append(etot)

    etots=np.array(etots)
    sort_dirs=np.array(dirs)[etots.argsort()]

    index=0.0
    for idir in sort_dirs:

        pm,etot=popu.get_final_properties(basepath+os.sep+idir)
        ea=np.array(pm['final_dmat']['euler_angles'])
        #print(idir,etot)
        #print(ea.shape)
        nangles=ea.shape[1]

        for j in range(nangles):
            theta = ea[:,j]
            dim=len(theta)
            radii = np.ones(dim)
            colors = np.arange(dim)
            width = 0.1*np.ones(dim)

            ax = fig.add_axes([float(j)/nangles, index/(len(dirs)+1), 1.0/nangles, 1.0/nangles], projection='polar')
            ax.yaxis.set_tick_params(labelsize=0)
            ax.xaxis.set_tick_params(labelsize=0)
            ax.get_xaxis().set_visible(False)
            ax.get_yaxis().set_visible(False)
            ax.spines['polar'].set_visible(True)
            bars = ax.bar(theta, radii, width=width, bottom=0.0)

            # Use custom colors and opacity
            for r, bar in zip(colors, bars):
                bar.set_facecolor(plt.cm.viridis(float(r)/nangles))
                bar.set_alpha(0.9)

        index+=1.0

    plt.savefig('OrbitalPolar.pdf')        
    plt.show()


def plot_status(basepath):

    dirs=[ x for x in os.listdir(basepath) 
           if os.path.isdir(basepath+os.sep+x) and os.path.isfile(basepath+os.sep+x+os.sep+'abinit.in') ]

    abi=AbinitInput(basepath+os.sep+'abinit.in')
    nstep = abi['nstep']

    print('Plotting Status...')
    X={}
    Y={}

    # Get the data:
    max_nruns=0
    for path in dirs:
        X[path]=np.array([])
        Y[path]={}
        Y[path]['nres2']=np.array([])
        Y[path]['etot']=np.array([])
        Y[path]['delta']=np.array([])

        for i in range(100):
            if os.path.isfile(basepath+os.sep+path+os.sep+'abinit_%02d.txt' % i):
                abo=AbinitOutput(basepath+os.sep+path+os.sep+'abinit_%02d.txt' % i)
                if not abo.is_finished:
                    print('This output is not finished')
                    continue
                if max_nruns < i:
                    max_nruns = i
                try:
                    energ = abo.get_energetics()
                except:
                    raise RuntimeError("Failed procesing: %s" % (basepath+os.sep+path+os.sep+'abinit_%02d.txt' % i))

                nres2 = energ['nres2'][-1]
                etot  = energ['etot'][-1]
                nscf  = len(energ['etot'])
                x = np.arange(nstep*i, nstep*i+ len(energ['nres2']))
                yetot = np.array(energ['etot'])
                ynres2 = np.array(energ['nres2'])
                ydelta = np.array(np.abs(energ['deltaEh']))
                X[path]=np.concatenate((X[path],x))
                Y[path]['etot']=np.concatenate((Y[path]['etot'],yetot))
                Y[path]['nres2']=np.concatenate((Y[path]['nres2'],ynres2))
                Y[path]['delta']=np.concatenate((Y[path]['delta'],ydelta))
                print("%s ETOT:%15.6f NRES2=%15.6e Num SCF=%3d" % (path,etot, nres2, nscf))

    # RESIDUAL
    plt.figure(figsize=(8,11))
    miny=1E-1
    maxy=1E-16    
    plt.subplots_adjust(left=0.1, bottom=0.05, right=0.99, top=0.99,
                    wspace=None, hspace=None)
    for path in dirs:
        if len(X)>0:
            plt.semilogy(X[path], Y[path]['nres2'], label=path[-4:], lw=0.1)
            if miny> min(Y[path]['nres2']):
                miny=min(Y[path]['nres2'])
            if maxy< max(Y[path]['nres2']):
                maxy=max(Y[path]['nres2'])

    plt.ylim(miny,maxy)

    for i in nstep*np.arange(max_nruns+1):
        plt.semilogy([i,i],[miny,maxy],'0.5', lw=0.1)

    plt.xlim(0,max([max(X[path]) for path in dirs]))
    plt.legend()
    plt.xlabel("SCF iteration")
    plt.ylabel("Density Residual$^2$")
    plt.savefig('Orbital_Residual.pdf')

    # ETOT
    miny=1E6
    maxy=-1E6
    avg=0
    minlen=100
    plt.figure(figsize=(8,11))
    for path in dirs:
        if len(X)>0:
            plt.plot(X[path], Y[path]['etot'], label=path[-4:])
            if len(Y[path]['etot'])<minlen:
                minlen=len(Y[path]['etot'])
                avg=np.average(Y[path]['etot'][-int(minlen/2):])
            if miny> min(Y[path]['etot'][-int(minlen/2):]):
                miny=min(Y[path]['etot'][-int(minlen/2):])
            if maxy< max(Y[path]['etot'][-int(minlen/2):]):
                maxy=max(Y[path]['etot'][-int(minlen/2):])

    plt.subplots_adjust(left=0.15, bottom=0.05, right=0.95, top=0.95,
                    wspace=None, hspace=None)

    newminy=miny-0.1*(maxy-miny)
    newmaxy=maxy+0.1*(maxy-miny)
    miny=newminy
    maxy=newmaxy

    plt.ylim(miny,maxy)
    for i in nstep*np.arange(max_nruns+1):
        plt.plot([i,i],[miny,maxy],'0.5',lw=0.1)

    plt.xlim(0,max([max(X[path]) for path in dirs]))
    plt.legend()
    plt.xlabel("SCF iteration")
    plt.ylabel("Energy")
    plt.savefig('Orbital_ETotal.pdf')

    # deltaEh
    plt.figure(figsize=(8,11))
    miny=1E-1
    for path in dirs:
        if len(X)>0:
            plt.semilogy(X[path], Y[path]['delta'], label=path[-4:], lw=0.1)
            if miny>min(Y[path]['delta']):
                miny=min(Y[path]['delta'])

    plt.subplots_adjust(left=0.1, bottom=0.05, right=0.99, top=0.99,
                    wspace=None, hspace=None)

    for i in nstep*np.arange(max_nruns+1):
        plt.semilogy([i,i],[miny,1E3],'0.5', lw=0.1)

    plt.ylim(miny,1E3)
    plt.xlim(0,max([max(X[path]) for path in dirs]))
    plt.legend()
    plt.xlabel("SCF iteration")
    plt.ylabel("Delta Energy")
    plt.savefig('Orbital_deltaEh.pdf')

def create_population(num_candidates):
    popu.random_population(num_candidates)

def prepare_folders(scrpath):
    for i in popu.members:
        popu.prepare_folder(i, workdir=scrpath, source_dir=scrpath)
    
if __name__=="__main__":

    parser = argparse.ArgumentParser(description='Orbital DFTU Evaluator and Analysis Tool')

    subparsers = parser.add_subparsers(help='commands', dest='subparser_name')

    # The create command
    create_parser = subparsers.add_parser('create', help='Create the database')
    
    create_parser.add_argument('-db_settings', type=str, 
                               help='Filename with PyChemiaDB settings to connect to server (JSON file)', 
                               required=False, default='dbsettings.json')

    # The populate command
    populate_parser = subparsers.add_parser('populate', help='Add candidates to the population (used for testing)')
    
    populate_parser.add_argument('-db_settings', type=str, 
                               help='Filename with PyChemiaDB settings to connect to server (JSON file)', 
                               required=False, default='dbsettings.json')

    populate_parser.add_argument('-size', type=int, 
                               help='Number of candidates to populate (default: 16)', 
                               required=False, default=16)
    populate_parser.add_argument('-basepath', type=str, 
                            help='Path where calculations are performed', 
                            required=False, default='.')

    # A run command
    run_parser = subparsers.add_parser('run', help='Run Evaluator')

    run_parser.add_argument('-db_settings', type=str, 
                            help='Filename with PyChemiaDB settings to connect to server (JSON file)', 
                            required=False, default='dbsettings.json')
    run_parser.add_argument('-pbs_settings', type=str, 
                            help='Filename with PBS settings for launching jobs', 
                            required=False, default='pbssettings.json')
    run_parser.add_argument('-basepath', type=str, 
                            help='Path where calculations are performed (default: current folder)', 
                            required=False, default='.')

    # The plot command
    plot_parser = subparsers.add_parser('plot', help='Generate several plots')

    plot_parser.add_argument('-db_settings', type=str, 
                             help='Filename with PyChemiaDB settings to connect to server (JSON file)', 
                             required=False, default='dbsettings.json')

    plot_parser.add_argument('-basepath', type=str, 
                     help='Path where calculations are performed', 
                     required=False, default='.')

    args = parser.parse_args()
    
    # Loading the settings to access the PyChemia Database
    if not os.path.isfile(args.db_settings):
        print("Could not read a PyChemiaDB settings file (JSON): %s" % args.db_settings)
        parser.print_help()
        sys.exit(1)
    rf=open(args.db_settings)
    db_settings=json.load(rf)
    rf.close()
    print("DB settings: %s" % db_settings )

    pcdb = get_database(db_settings)

    if args.subparser_name == 'run':

        # Loading the settings to access to submit jobs with PBS
        if not os.path.isfile(args.pbs_settings):
            print("Could not read a PBS settings file (JSON): %s" % args.pbs_settings)
            parser.print_help()
            sys.exit(1)
        rf=open(args.pbs_settings)
        pbs_settings=json.load(rf)
        rf.close()
        print("PBS settings: %s" % pbs_settings )
    
    if args.subparser_name in ['run', 'plot', 'populate']:
        
        if not os.path.isdir(args.basepath) or not os.path.isfile(args.basepath+'/abinit.in'):
            print('ERROR: Wrong basepath %s, directory must exist and contain a abinit.in file' % args.basepath)
            parser.print_help()
            sys.exit(1)

        popu = OrbitalDFTU(pcdb, args.basepath+'/abinit.in')

    if args.subparser_name == 'populate':
        print(popu)
        popu.random_population(args.size)
        print("Total number of candidates: %d" % len(popu))

    if args.subparser_name == 'run':
        popu.evaluator(pbs_settings=pbs_settings, basedir=args.basepath)

    if args.subparser_name == 'plot':        
        check_status(args.basepath)
        plot_status(args.basepath)
        plot_polar(popu, args.basepath)
