import pychemia
import tempfile
from .samples import CaTiO3


def test_xyz():
    """
    Test (pychemia.io.xyz)                                      :
    """
    st1 = CaTiO3()
    st1.set_periodicity(False)
    file = tempfile.NamedTemporaryFile()
    pychemia.io.xyz.save(st1, file.name)
    st2 = pychemia.io.xyz.load(file.name)
    assert st1 == st2


def test_ascii():
    """
    Test (pychemia.io.ascii)                                    :
    """
    st1 = CaTiO3()
    file = tempfile.NamedTemporaryFile()
    pychemia.io.ascii.save(st1, file.name)
    st2 = pychemia.io.ascii.load(file.name)
    return st1, st2
