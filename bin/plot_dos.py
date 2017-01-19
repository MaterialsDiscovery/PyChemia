#!/usr/bin/env python

import os
import sys
import pychemia


def helper():
    print(""" Plot the density of states from a given list of files
   Use:

       plot_dos.py [--figname 'DensityOfStates.pdf' ]
                   [--set_energy_range min max ]
                   [--set_figsize figwidth figheight]
                   [--help ]
                   file1.dat file2.dat ...
   """)


if __name__ == '__main__':

    # Script starts from here
    if len(sys.argv) < 2:
        helper()
        sys.exit(1)

    figname = 'DensityOfStates.pdf'
    minenergy = None
    maxenergy = None
    figheight = None
    figwidth = None
    filelist = []
    doslist = []

    for i in range(1, len(sys.argv)):
        if sys.argv[i].startswith('--'):
            option = sys.argv[i][2:]
            # fetch sys.argv[1] but without the first two characters
            if option == 'figname':
                figname = sys.argv[i + 1]
            elif option == 'help':
                helper()
                sys.exit()
            elif option == 'set_energy_range':
                minenergy = float(sys.argv[i + 1])
                maxenergy = float(sys.argv[i + 2])
            elif option == 'set_figsize':
                figwidth = float(sys.argv[i + 1])
                figheight = float(sys.argv[i + 2])
            else:
                print('Unknown option. --' + option)
        elif os.path.isfile(sys.argv[i]) and sys.argv[i][-4:] == '.dat':
            filelist.append(sys.argv[i])

    print(figname)
    for i in filelist:
        a = pychemia.visual.DensityOfStates().read(i)
        doslist.append(a)

    if len(doslist) == 0:
        helper()
        sys.exit(1)
    elif len(doslist) == 1:
        fig, ax = pychemia.visual.plot_one_dos(doslist[0], ax=None, horizontal=True, figwidth=figwidth,
                                               figheight=figheight)
        # fig.savefig(figname)
    else:
        fig, ax = pychemia.visual.plot_many_dos(doslist, minenergy=minenergy, maxenergy=maxenergy, figwidth=figwidth,
                                                figheight=figheight)
        # fig.savefig(figname)
