__author__ = 'Guillermo Avendano Franco'


def test_symmetry():
    """
    Tests from doctests for symmetry    :
    """
    import doctest
    import pychemia.symm

    dt = doctest.testmod(pychemia.symm.symmetry, verbose=True)
    assert dt.failed == 0
