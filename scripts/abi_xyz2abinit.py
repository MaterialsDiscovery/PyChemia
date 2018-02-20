#!/usr/bin/env python
import os
import sys
import pychemia.code.abinit


if __name__ == '__main__':

    # Script starts from here
    if len(sys.argv) < 2:
        print('No file specified.')
        sys.exit()

    filename = sys.argv[1]
    assert (os.path.isfile(filename))
    iv = pychemia.code.abinit.xyz2input(filename)
    av = pychemia.code.abinit.AbinitInput()
    av.variables = iv.variables
    av.write(filename + ".in")
