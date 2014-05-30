__author__ = 'Guillermo Avendano-Franco'


def test_periodic():
    """
    Tests from doctests for periodic    :
    """
    import doctest
    import pychemia.utils.periodic

    doctest.testmod(pychemia.utils.periodic, verbose=True)


def test_mathematics():
    """
    Tests from doctests for mathematics :
    """
    import doctest
    import pychemia.utils.mathematics

    doctest.testmod(pychemia.utils.mathematics, verbose=True)
