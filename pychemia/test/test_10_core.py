import pychemia
import numpy as np


def test_composition():
    """
    Test class pychemia.Composition     :
    """
    comp = pychemia.Composition('YBa2Cu3O7')
    assert (comp.species == sorted([u'Y', u'Cu', u'Ba', u'O']))
    assert (
    sorted(comp.symbols) == sorted([u'Y', u'Cu', u'Cu', u'Cu', u'Ba', u'Ba', u'O', u'O', u'O', u'O', u'O', u'O', u'O']))
    assert (comp.formula == 'Ba2Cu3O7Y')
    assert (abs(comp.covalent_volume() - 285.185) < 1E-5)
    assert (comp.natom == 13)
    assert (comp.species_hex() == 942087432)
    assert (comp.species_bin() == 72058144330612992)
    assert (comp.sorted_formula(sortby='hill') == 'Ba2Cu3O7Y')

    comp = pychemia.Composition('Na2Cl2')
    assert (sorted(comp.symbols) == sorted([u'Na', u'Na', u'Cl', u'Cl']))
    assert (pychemia.utils.periodic.valence(comp.symbols) == [1, 1, 7, 7])


def test_structure():
    """
    Test class pychemia.Structure       :
    """
    a = 4.05
    b = a / 2
    fcc = pychemia.Structure(symbols=['Au'], cell=[[0, b, b], [b, 0, b], [b, b, 0]], periodicity=True)
    assert (fcc.natom == 1)
    fcc_copy = fcc.copy()
    fcc_copy.canonical_form()
    assert (abs(fcc.volume - fcc_copy.volume) < 1E-13)
    assert (np.linalg.norm(fcc_copy.lattice.angles - fcc.lattice.angles) < 1E-10)
    spc = fcc.supercell((3, 3, 3))
    assert (spc.natom == 27)
