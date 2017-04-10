#!/usr/bin/env python

from builtins import input

import os
import sys
import argparse
import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import pychemia

def compare_params(path):

    abinitout = get_final_abinit_out(path)

    if not os.path.isfile(path+os.sep+'abinit.in'):
        print('ERROR: No abinit.in found at %s' % path)
        return

    # For making easier to see the values
    np.set_printoptions(linewidth=200, suppress=True)

    # Reading the INPUT
    abi=pychemia.code.abinit.InputVariables(path+os.sep+'abinit.in')
    idmatpawu=np.array(abi['dmatpawu']).reshape(-1,5,5)
    iparams=pychemia.population.orbitaldftu.dmatpawu2params(idmatpawu,5)

    # Reading the OUTPUT
    dmatpawu=pychemia.population.orbitaldftu.get_final_dmatpawu(abinitout)
    odmatpawu=np.array(dmatpawu).reshape(-1,5,5)
    oparams=pychemia.population.orbitaldftu.dmatpawu2params(odmatpawu,5)

    #print('DMATPAWU')
    #print('input')
    #print(idmatpawu)
    #print('output')
    #print(odmatpawu)

    print('PARAMETRIC REPRESENTATION found at %s' % abinitout)
    for i in sorted(list(iparams.keys())):
        print(i)
        print('input')
        print(iparams[i])
        print('output')
        print(i)
        print(oparams[i])

    abo=pychemia.code.abinit.AbinitOutput(abinitout)
    if not abo.is_finished:
        print('This output is not finished')
    try:
        nres2=abo.get_energetics()['nres2'][-1]
        etot=abo.get_energetics()['etot'][-1]
        nscf=len(abo.get_energetics()['etot'])
        print("%30s ETOT: %15.6f NRES2: %15.6e NumSCF: %3d" % (path,etot, nres2, nscf))
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
        path=basepath+os.sep+i
        abinitout = get_final_abinit_out(path)
        if abinitout is None:
            continue

        abo=pychemia.code.abinit.AbinitOutput(abinitout)
        if not abo.is_finished:
            continue
        try:
            nres2=abo.get_energetics()['nres2'][-1]
            etot=abo.get_energetics()['etot'][-1]
            nscf=len(abo.get_energetics()['etot'])
            print("%-40s %15.6f %15.6e %4d" % (abinitout, etot, nres2, nscf))
        except:
            print("ERROR: Could not get final energetics from %s" % (abinitout))

def plot_status(basepath,dirs, nstep=50):
    print('Plotting Status...')
    X={}
    Y={}

    # Get the data:
    for path in dirs:
        X[path]=np.array([])
        Y[path]={}
        Y[path]['nres2']=np.array([])
        Y[path]['etot']=np.array([])
        Y[path]['delta']=np.array([])

        for i in range(10):
            if os.path.isfile(basepath+os.sep+path+os.sep+'abinit_%d.out' % i):
                abo=pychemia.code.abinit.AbinitOutput(basepath+os.sep+path+os.sep+'abinit_%d.out' % i)
                if not abo.is_finished:
                    print('This output is not finished')
                try:
                    energ = abo.get_energetics()
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
                except:
                    print("%s Failed processing output" % path)

    # RESIDUAL
    plt.figure(figsize=(8,11))
    plt.subplots_adjust(left=0.1, bottom=0.05, right=0.99, top=0.99,
                    wspace=None, hspace=None)
    for path in dirs:
        if len(X)>0:
            plt.semilogy(X[path], Y[path]['nres2'], label=path[-4:])

    for i in nstep*np.arange(10):
        plt.semilogy([i,i],[1E-19,1E5],'0.5')

    plt.ylim(1E-19,1E5)
    plt.xlim(0,max([max(X[path]) for path in dirs]))
    plt.legend()
    plt.xlabel("SCF iteration")
    plt.ylabel("Density Residual$^2$")
    plt.savefig('Residual.pdf')

    # ETOT
    plt.figure(figsize=(8,11))
    for path in dirs:
        if len(X)>0:
            plt.plot(X[path], Y[path]['etot'], label=path[-4:])

    plt.subplots_adjust(left=0.1, bottom=0.05, right=0.99, top=0.99,
                    wspace=None, hspace=None)

    for i in nstep*np.arange(10):
        plt.plot([i,i],[-992,-990],'0.5')

    plt.ylim(-992,-990)
    plt.xlim(0,max([max(X[path]) for path in dirs]))
    plt.legend()
    plt.xlabel("SCF iteration")
    plt.ylabel("Energy")
    plt.savefig('ETotal.pdf')

    # deltaEh
    plt.figure(figsize=(8,11))
    for path in dirs:
        if len(X)>0:
            plt.semilogy(X[path], Y[path]['delta'], label=path[-4:])

    plt.subplots_adjust(left=0.1, bottom=0.05, right=0.99, top=0.99,
                    wspace=None, hspace=None)

    for i in nstep*np.arange(10):
        plt.semilogy([i,i],[1E-13,1E3],'0.5')

    plt.ylim(1E-13,1E3)
    plt.xlim(0,max([max(X[path]) for path in dirs]))
    plt.legend()
    plt.xlabel("SCF iteration")
    plt.ylabel("Delta Energy")
    plt.savefig('deltaEh.pdf')


def plot_polar(basepath,dirs):
    print('Plotting Euler Angles...')

    fig=plt.figure(figsize=(10,1.2*len(dirs)))
    plt.subplots_adjust(left=0.0, bottom=0.0, right=0.99, top=0.99, wspace=None, hspace=None)
    etots=[]

    for idir in dirs:
        pm, etot=get_final_properties(basepath+os.sep+idir)
        etots.append(etot)

    etots=np.array(etots)
    sort_dirs=np.array(dirs)[etots.argsort()]

    index=0.0
    for idir in sort_dirs:

        pm,etot=get_final_properties(basepath+os.sep+idir)
        ea=np.array(pm['final_dmat']['euler_angles'])
        print(idir,etot)

        for j in range(10):
            theta = ea[:,j]
            radii = np.ones(4)
            colors = np.arange(4)
            width = 0.1*np.ones(4)

            ax = fig.add_axes([float(j)/10, index/(len(dirs)+1), 0.1, 0.1], projection='polar')
            ax.yaxis.set_tick_params(labelsize=0)
            ax.xaxis.set_tick_params(labelsize=0)
            ax.get_xaxis().set_visible(False)
            ax.get_yaxis().set_visible(False)
            ax.spines['polar'].set_visible(False)
            bars = ax.bar(theta, radii, width=width, bottom=0.0)

            # Use custom colors and opacity
            for r, bar in zip(colors, bars):
                bar.set_facecolor(plt.cm.viridis(r / 10.))
                bar.set_alpha(0.9)

        index+=1.0

    plt.savefig('OrbitalPolar.pdf')        
    plt.show()
    
if __name__=="__main__":

    parser = argparse.ArgumentParser(description='Orbital DFTU Evaluator and Analysis Tool')

    parser.add_argument('-host', type=str, help='Hostname or IP of the Mongo Server', required=True)
    parser.add_argument('-name', type=str, help='Name of Database', required=True)
    parser.add_argument('-dbuser', type=str, help='Username on Database', required=True)
    parser.add_argument('-quser', type=str, help='Username on PBS Queue ', required=True)
    parser.add_argument('-queue', type=str, help='Queue', required=True)
    parser.add_argument('-hours', type=int, help='Number of hours to request', required=True)
    parser.add_argument('-features', type=str, help='Features for queue', required=False, default=None)
    parser.add_argument('-basepath', type=str, help='Path where calculations are performed', required=False, default='.')

    args = parser.parse_args()
    basepath = args.basepath    

    passwd = input('Password:')

    if not os.path.isdir(args.basepath) or not os.path.isfile(basepath+'/abinit.in'):
        print('ERROR: Wrong basepath %s' % basepath)
        parser.print_help()
        sys.exit(1)

    db_settings={'name': args.name, 'host': args.host, 'user':args.dbuser, 'passwd': passwd}
    pcdb = pychemia.db.get_database(db_settings)
    popu=pychemia.population.OrbitalDFTU(pcdb, basepath+'/abinit.in')
    popu.evaluator(username=args.quser, basedir=basepath, queue=args.queue, walltime=[args.hours,0,0], ppn=6, features=args.features)
