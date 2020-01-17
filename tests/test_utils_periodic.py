from pychemia.utils.periodic import *


def test_utils():
    """
    Test (pychemia.utils.periodic)                              :
    """
    assert covalent_radius(1) == 0.31
    assert covalent_radius([1.1, 2]) == [0.31, 0.28]
    assert len(covalent_radius()) == 118
    assert atomic_symbol(1) == 'H'
    assert atomic_symbol([1, 2.1]) == ['H', 'He']
    assert len(atomic_symbol()) == 118


def test_atomicnumber():
    """
    Test (pychemia.utils.periodic) [atomic_number]              :
    """
    assert (atomic_number('H') == 1)
    assert (atomic_number(['H', 'He']) == [1, 2])


def test_mass():
    """
    Test (pychemia.utils.periodic) [mass]                       :
    """
    assert (mass(['H', 'He']) == [1.00794, 4.002602])
    assert (mass([1, 2]) == [1.00794, 4.002602])
    assert (mass(1) == 1.00794)


def test_symbol():
    """
    Test (pychemia.utils.periodic) [atomic_symbol]              :
    """
    assert (atomic_symbol(1) == 'H')
    assert (atomic_symbol([1, 2]) == ['H', 'He'])


def test_covalent_radius():
    """
    Test (pychemia.utils.periodic) [covalent_radius]            :
    """
    assert (covalent_radius('H') == 0.31)
    assert (covalent_radius([1, 2]) == [0.31, 0.28])
    assert (covalent_radius(['H', 'He']) == [0.31, 0.28])


def test_valence():
    """
    Test (pychemia.utils.periodic) [valence]                    :
    """
    assert (valence(1) == 1)
    assert (valence([1, 2]) == [1, 0])
