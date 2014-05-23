#!/usr/bin/env python

# This script download several sets of pseudopotentials
# form the Website of ABINIT
# The pseudopotentials will be located inside the directory
# '.abinit' of the user home directory, just aside of the
# tarballs.
# AE are the Hirshfeld All Electron (AE) densities

import os

import pychemia.dft.codes.abinit as pa


if __name__ == '__main__':

    home = os.environ['HOME']
    basedir = home + "/.abinit"

    if not os.path.isdir(basedir):
        os.mkdir(basedir)

    for exchange in ['LDA', 'GGA', 'AE']:

        lista = []
        if exchange == 'LDA':
            lista = ['FHI', 'TM']
        elif exchange == 'GGA':
            lista = ['FHI']
        elif exchange == 'AE':
            lista = ['DEN']

        for kind in lista:
            pa.get_all_psps(basedir, exchange, kind)
