#!/usr/bin/env python

import os
import sys
import pychemia.code.octopus


def helper():
    print(""" Set a variable in octopus
   Use:

       octopus_setvar.py --filename 'Octopus_Input_File' [ --set varname value ]...  [--del varname]...
   """)


if __name__ == '__main__':

    # Script starts from here
    if len(sys.argv) < 2:
        helper()
        sys.exit(1)

    filename = ''
    toset = {}
    todel = []
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
            elif option == 'filename':
                filename = sys.argv[i + 1]
            elif option == 'set':
                toset[sys.argv[i + 1]] = sys.argv[i + 2]
            elif option == 'del':
                todel.append(sys.argv[i + 1])
            else:
                print('Unknown option. --' + option)

    # Set the variables
    if os.path.isfile(filename):
        data = pychemia.code.octopus.AbinitInput(filename)
        for i in toset.keys():
            data.variables[i] = toset[i]
        for i in todel:
            if i in data.variables.keys():
                data.variables.pop(i)
        data.write(filename)
    else:
        print('ERROR: no filename', filename)
