__author__ = 'Guillermo'

import os
import pychemia


def test_incar():
    """
    Test VASP INCAR parsing and writing :
    """
    print os.getcwd()
    vi = pychemia.code.vasp.read_incar('pychemia/test/data/INCAR_test')
    rf = open('pychemia/test/data/INCAR_chk')
    assert (str(vi) == rf.read())

if __name__ == '__main__':
    test_incar()
