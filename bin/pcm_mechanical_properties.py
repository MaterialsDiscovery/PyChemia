import getopt
import logging
import os
import sys

import numpy as np

from pychemia.code.vasp.task.elastic import mechanical_properties, elastic_moduli

logging.basicConfig(level=logging.DEBUG)


def usage(name):
    print("""
NAME
    %s

DESCRIPTION
    Compute the mechanical properties from a given 'OUTCAR' file

OPTIONS

    --help, -h
        Return information on the options and use of this script

    --outcar, -o <string> (Default: 'OUTCAR')
        OUTCAR file from which the mechanical properties will be computed

""" % os.path.basename(name))


def main(argv):
    try:
        opts, args = getopt.getopt(argv[1:], "ho:", ["help", "outcar="])
    except getopt.GetoptError:
        usage(argv[0])
        sys.exit(2)

    # Default Values
    outcar = 'OUTCAR'

    for opt, arg in opts:

        if opt in ("-h", "--help"):
            usage(argv[0])
            sys.exit()
        elif opt in ("-o", "--outcar"):
            outcar = arg

    tem = np.array(elastic_moduli(outcar)['total_elastic_moduli'])
    np.set_printoptions(precision=3)
    print(' Total Elastic Moduli')
    print(' --------------------\n')
    for irow in tem:
        for jelem in irow:
            if jelem == 0.0 or jelem == -0.0:
                sys.stdout.write(' %8d    ' % 0)
            else:
                sys.stdout.write(' %12.3f' % jelem)
        sys.stdout.write('\n')

    print('\n')
    mech = mechanical_properties(tem)
    print(" %20s %5s %12s %12s %12s" % ('Property'.ljust(20), 'Units', 'Voigt', 'Reuss', 'Average'))
    print(' ' + 65 * "-")
    for i in sorted(mech):
        print(" %20s %5s %12.3f %12.3f %12.3f" % (i.ljust(20), mech[i]['units'], mech[i]['Voigt'], mech[i]['Reuss'],
                                                  0.5 * (mech[i]['Voigt'] + mech[i]['Reuss'])))
    print(' ' + 65 * "-")


if __name__ == "__main__":
    main(sys.argv)
