__author__ = 'Guillermo Avendano Franco'


def test_symmetry():
    """
    Tests from doctests for symmetry    :
    """
    import doctest
    import pychemia
    if pychemia.HAS_SPGLIB:
        import pychemia.symm
        dt = doctest.testmod(pychemia.symm.symmetry, verbose=True)
        assert dt.failed == 0
