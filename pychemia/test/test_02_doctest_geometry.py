__author__ = 'Guillermo'


def test():
    """
    Tests from doctests for structure   :
    """
    import doctest
    import pychemia.geometry.structure

    doctest.testmod(pychemia.geometry.structure, verbose=True)
