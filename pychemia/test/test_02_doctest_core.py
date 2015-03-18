__author__ = 'Guillermo'


def test_structure():
    """
    Tests from doctests for structure   :
    """
    import doctest
    import pychemia.core.structure

    doctest.testmod(pychemia.core.structure, verbose=True)


def test_lattice():
    """
    Tests from doctests for lattice     :
    """
    import doctest
    import pychemia.core.lattice

    doctest.testmod(pychemia.core.lattice, verbose=True)
