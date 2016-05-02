#!/usr/bin/env python
import os
import sys
from pychemia import HAS_SCIPY

if HAS_SCIPY:
    import pychemia.code.abinit
else:
    raise ImportError('scipy could not be found')

if __name__ == '__main__':

    # Script starts from here
    if len(sys.argv) < 2:
        print('No file specified.')
        sys.exit()

    filename = sys.argv[1]
    assert (os.path.isfile(filename))
    iv = pychemia.code.abinit.xyz2input(filename)
    av = pychemia.code.abinit.InputVariables()
    av.variables = iv.variables
    av.write(filename + ".in")
