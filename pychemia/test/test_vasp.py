import os
import pychemia


def test_incar():
    """
    Test VASP INCAR parsing and writing :
    """
    print os.getcwd()
    vi = pychemia.code.vasp.read_incar('pychemia/test/data/vasp_05/INCAR_test')
    rf = open('pychemia/test/data/vasp_05/INCAR_chk')
    assert (str(vi) == rf.read())

if __name__ == '__main__':
    test_incar()


def test_bad_outcar():

    vo = pychemia.code.vasp.VaspOutput('pychemia/test/data/VASP_04/OUTCAR')
    assert not vo.is_finished