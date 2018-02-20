#!/usr/bin/env python

import os
import sys
import numpy as np
from pychemia.utils.constants import bohr_angstrom
import pychemia.code.abinit


def input2xyz(abivar, basename, datasets):
    """
    Take an input object 'av' and write the
    corresponding XYZ file in the given
    'basename'
    """

    for idts in datasets:
        struct = abivar.get_structure(idts)
        natom = struct.natom
        if idts == '':
            filep = basename + '.xyz'
        else:
            filep = basename + '_DS' + str(idts) + '.xyz'
        wf = open(filep, 'w')
        wf.write(str(natom) + '\n\n')
        for iatom in range(natom):
            wf.write('%s %16.8f %16.8f %16.8f\n' % (struct.symbols[iatom],
                                                    struct.positions[iatom, 0],
                                                    struct.positions[iatom, 1],
                                                    struct.positions[iatom, 2]))
        wf.close()


def write_one(time, nchist, wf, natom, struct):
    xcart = np.array(nchist['xcart'][time])
    atom_pos = bohr_angstrom * xcart
    wf.write(str(natom) + '\n\n')
    for iatom in range(natom):
        wf.write('%s %16.8f %16.8f %16.8f\n' % (struct['atom_names'][iatom], atom_pos[iatom, 0], atom_pos[iatom, 1],
                                                atom_pos[iatom, 2]))


def abihist2xyz(abivar, basename, datasets, time='all'):
    """
    Write XYZ files for history
    """
    for idts in datasets:

        # Reading the history
        if idts == '':
            history = abifile.files['tmpout'] + "_HIST"
        else:
            history = abifile.files['tmpout'] + "_DS" + str(idts) + "_HIST"
        if os.path.isfile(history):
            print('Reading ', history)
            nchist = pychemia.code.abinit.netcdf2dict(history)

            # Setting the output file
            if time == 'all':
                if idts == '':
                    filep = basename + '_HIST.xyz'
                else:
                    filep = basename + '_DS' + str(idts) + '_HIST.xyz'
            else:
                if idts == '':
                    filep = basename + '_HIST_itime' + str(time) + '.xyz'
                else:
                    filep = basename + '_DS' + str(idts) + '_HIST_itime' + str(time) + '.xyz'

            print('Writing ', filep)
            wf = open(filep, 'w')

            # Getting the atomic structure
            struct = abivar.get_struct(idts)
            ntime = len(nchist['mdtime'])
            natom = struct['natom']

            # Write the xyz section
            if time == 'all':
                sys.stdout.write('[')
                for it in range(ntime):
                    sys.stdout.write('*')
                    write_one(it, nchist, wf, natom, struct)
                sys.stdout.write(']\n')
            else:
                write_one(time, nchist, wf, natom, struct)

            # Close the file
            wf.close()


def helper():
    print('''\
   Extract atomic coordinates from abinit files
   This program can deal with the following
   files:
           abinit.files ('.files': abinit format)
           abinit.in    ('.in': input for abinit)
           *_OUT.nc     ('_OUT.nc': output in NetCDF file format)

   Options include:
     --version : Prints the version number
     --help    : Display this help
     --dtset # : Extract from a given dataset
     --itime # : Extract only a given iteration (--hist mandatory)
     --hist    : Extract from HIST file
     --input   : Extract from input file (default)
     --output  : Extract from output file
     ''')


def get_input_dts(filep, dataset):
    abivar = pychemia.code.abinit.AbinitInput(filep)
    datasets = None
    if dataset == 0:
        datasets = abivar.get_dtsets_keys()
    elif not str(dataset) in abivar.get_dtsets_keys():
        print('ERROR No such dataset index:', dataset)
        print('Valid values are:', abivar.get_dtsets_keys())
        exit(1)
    else:
        datasets = [dataset]
    return abivar, datasets


if __name__ == '__main__':

    # Script starts from here
    if len(sys.argv) < 2:
        helper()
        sys.exit(1)

    filename = ''
    dtset = 0
    itime = 'all'
    has_input = True
    has_output = False
    hist = False
    for i in range(1, len(sys.argv)):
        if sys.argv[i].startswith('--'):
            option = sys.argv[i][2:]
            # fetch sys.argv[1] but without the first two characters
            if option == 'version':
                print('Version 1.0')
                sys.exit()
            elif option == 'help':
                helper()
                sys.exit()
            elif option == 'dtset':
                dtset = int(sys.argv[i + 1])
            elif option == 'itime':
                itime = int(sys.argv[i + 1])
            elif option == 'input':
                has_input = True
                has_output = False
            elif option == 'output':
                has_output = True
                has_input = False
            elif option == 'hist':
                has_output = False
                has_input = False
                hist = True
            else:
                print('Unknown option. --' + option)

    filename = sys.argv[-1]
    if not os.path.isfile(filename):
        print('No valid filename')
        helper()
        sys.exit(1)

    av = None
    dts = None
    if filename != '':
        if filename[-6:] == '.files':
            abifile = pychemia.code.abinit.AbiFiles(filename)
            if has_input or hist:
                (av, dts) = get_input_dts(abifile.get_input_filename(), dtset)
            elif has_output:
                (av, dts) = get_input_dts(abifile.get_out_filename(), dtset)
            xyz = filename[:-6]
            if not hist:
                input2xyz(av, xyz, dts)
            else:
                abihist2xyz(av, xyz, dts, itime)
        elif filename[-3:] == '.in':
            (av, dts) = get_input_dts(filename, dtset)
            xyz = filename[:-3]
            input2xyz(av, xyz, dts)
        elif filename[-7:] == '_OUT.nc':
            (av, dts) = get_input_dts(filename, dtset)
            xyz = filename[:-7]
            input2xyz(av, xyz, dts)
        else:
            print('Filename with unknown extension')
    else:
        print('No file specified.')
        helper()
        sys.exit(1)
