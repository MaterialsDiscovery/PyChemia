import doctest


def test_periodic():
    """
    Tests from doctests for periodic    :
    """
    import pychemia.utils.periodic
    dt = doctest.testmod(pychemia.utils.periodic, verbose=True)
    assert dt.failed == 0


def test_mathematics():
    """
    Tests from doctests for mathematics :
    """

    import pychemia.utils.mathematics
    dt = doctest.testmod(pychemia.utils.mathematics, verbose=True)
    assert dt.failed == 0
