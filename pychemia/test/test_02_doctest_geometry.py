__author__ = 'Guillermo'


def test():
    """
    Tests from doctests for structure   :
    """
    import doctest
    import pychemia.core.structure

    doctest.testmod(pychemia.core.structure, verbose=True)
