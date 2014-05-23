__author__ = 'viviane'


def test():
    """
    Tests from doctests for periodic    :
    """
    import doctest
    import pychemia.utils.periodic

    doctest.testmod(pychemia.utils.periodic, verbose=True)
