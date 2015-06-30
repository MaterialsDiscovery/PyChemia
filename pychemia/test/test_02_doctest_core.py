__author__ = 'Guillermo'

import doctest


def test_structure():
    """
    Tests from doctests for structure   :
    """
    import pychemia.core.structure
    dt = doctest.testmod(pychemia.core.structure, verbose=True)
    assert dt.failed == 0


def test_lattice():
    """
    Tests from doctests for lattice     :
    """
    import pychemia.core.lattice
    dt = doctest.testmod(pychemia.core.lattice, verbose=True)
    assert dt.failed == 0


def test_composition():
    """
    Tests from doctests for composition :
    """
    import pychemia.core.composition
    dt = doctest.testmod(pychemia.core.composition, verbose=True)
    assert dt.failed == 0
