import sys

__author__ = 'Guillermo Avendano Franco'

import os
import numpy as np
import matplotlib.pyplot as plt

class DensityOfStates():

    def __init__(self, table = None, title = None):
        self._dos = None
        self._min_energy = None
        self._max_energy = None
        self._max_dos = None
        self.title = title
        self.ncols = 1

        if table is not None:
            self._dos = np.array(table)
            self.ncols= self._dos.shape[1]-1
            self._min_energy = min(table[:,0])
            self._max_energy = max(table[:,0])
            self._max_dos = max(table[:,1])

    @staticmethod
    def read(filename, title = None):

        table = np.loadtxt(filename)
        if title is None:
            root, ext = os.path.splitext(os.path.basename(filename))
            name = root
        else:
            name=title

        jump=0
        jumplist=[0]
        nsets=1
        for iline in range(len(table)-1):
            if table[iline,0]!=table[iline-jump,0]:
                print iline, jump, table[iline,0], table[iline-jump,0]
                raise ValueError("No consistency on energy values")
            if table[iline+1, 0] < table[iline, 0]:
                jump=iline+1
                jumplist.append(jump)
                nsets+=1

        if nsets>1:
            jump1=jumplist[1]-jumplist[0]
            for i in range(1,nsets-1):
                if jumplist[i+1]-jumplist[i]!=jump1:
                    raise ValueError("No equal jumps")
            table2=np.zeros((jump1,nsets+1))
            table2[:,0]=table[:jump1,0]
            for i in range(nsets):
                table2[:,i+1]=table[i*jump1:(i+1)*jump1,1]
                assert(np.all(table[i*jump1:(i+1)*jump1,0]==table[:jump1,0]))
        else:
            table2=table

        dos = DensityOfStates(table= table2, title= name)
        return dos

    @property
    def energies(self):
        return self._dos[:,0]

    @property
    def values(self):
        if self.ncols>1:
            return self._dos[:,range(1,self.ncols+1)]
        else:
            return self._dos[:,1]


def plot_one_dos(dosobj, ax=None, horizontal=True):

    if ax is None:
        fig = plt.figure()
        ax = plt.subplot(111)

    X = dosobj.energies

    if dosobj.ncols>1:
        for i in range(dosobj.ncols):
            Y = dosobj.values[:,i]

            if horizontal:
                ax.plot(X, Y)
            else:
                ax.plot(Y, X)

    else:
        Y = dosobj.values
        if horizontal:
            ax.plot(X, Y)
        else:
            ax.plot(Y, X)


def plot_many_dos(doslist, minenergy = None, maxenergy = None):
    ndos = len(doslist)
    if minenergy is None:
        minenergy = min([ min(x.energies) for x in doslist])
    if maxenergy is None:
        maxenergy = max([ max(x.energies) for x in doslist])
    minval = sys.float_info.max
    maxval = sys.float_info.min
    for idos in doslist:
        for i in idos._dos:
            if i[0]>minenergy and i[0]<maxenergy:
                for icol in range(idos.ncols):
                    if i[icol+1]>maxval:
                        maxval = i[icol+1]
                    if i[icol+1]<minval:
                        minval = i[icol+1]
    print minval, maxval

    fig, ax = plt.subplots(nrows=1, ncols=ndos, sharex=False, sharey=True, squeeze=True)
    fig.set_figheight(12)
    fig.set_figwidth(16)
    plt.subplots_adjust(left=0.05, bottom=0.05, right=0.95, top=0.95, wspace=0, hspace=0)
    for i in range(ndos):
        plot_one_dos(doslist[i], ax[i], horizontal= False)
        ax[i].set_xlim(1.1*minval,1.1*maxval)
        ax[i].set_ylim(minenergy,maxenergy)
        ax[i].set_xlabel(doslist[i].title)
    ax[0].set_ylabel('Energy')
    plt.savefig('dos.pdf')
