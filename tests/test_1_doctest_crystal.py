import doctest


def test_symmetry():
    """
    DocTests (pychemia.crystal.symmetry)                         :
    """
    import pychemia
    if pychemia.HAS_SPGLIB:
        import pychemia.crystal
        dt = doctest.testmod(pychemia.crystal.symmetry, verbose=True)
        assert dt.failed == 0


def test_lattice():
    """
    DocTests (pychemia.crystal.lattice)                          :
    """
    import pychemia.crystal.lattice
    dt = doctest.testmod(pychemia.crystal.lattice, verbose=True, optionflags=doctest.NORMALIZE_WHITESPACE)
    assert dt.failed == 0
