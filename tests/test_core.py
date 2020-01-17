import pychemia
import numpy as np
import unittest


class TestCore(unittest.TestCase):
    """
    Test (pychemia.core)                            :
    """

    def test_composition_1(self):
        """
        Test (pychemia.core.composition)                            :
        """
        comp = pychemia.Composition('YBa2Cu3O7')
        self.assertEqual(sorted(comp.species), [u'Ba', u'Cu', u'O', u'Y'])
        self.assertEqual(sorted(comp.symbols),
                         [u'Ba', u'Ba', u'Cu', u'Cu', u'Cu', u'O', u'O', u'O', u'O', u'O', u'O', u'O', u'Y'])
        self.assertEqual(comp.formula, 'Ba2Cu3O7Y')
        self.assertTrue(abs(comp.covalent_volume() - 285.185) < 1E-5)
        self.assertEqual(comp.natom, 13)
        self.assertEqual(comp.sorted_formula(sortby='hill'), 'Ba2Cu3O7Y')

    def test_composition_2(self):
        """
        Test (pychemia.core.composition)                            :
        """
        comp = pychemia.Composition('Na2Cl2')
        self.assertEqual(comp.symbols, [u'Cl', u'Cl', u'Na', u'Na'])
        self.assertEqual(pychemia.utils.periodic.valence(comp.symbols), [7, 7, 1, 1])

    def test_structure(self):
        """
        Test (pychemia.core.structure)                              :
        """
        a = 4.05
        b = a / 2
        fcc = pychemia.Structure(symbols=['Au'], cell=[[0, b, b], [b, 0, b], [b, b, 0]], periodicity=True)
        self.assertEqual(fcc.natom, 1)
        fcc_copy = fcc.copy()
        fcc_copy.canonical_form()
        self.assertTrue(abs(fcc.volume - fcc_copy.volume) < 1E-13)
        self.assertTrue(np.linalg.norm(fcc_copy.lattice.angles - fcc.lattice.angles) < 1E-10)
        spc = fcc.supercell((3, 3, 3))
        self.assertEqual(spc.natom, 27)

    def test_from_file_1(self):
        """
        Test (pychemia.core.from_file)                              :
        """
        filename = 'tests/data/vasp_07/POSCAR_new'
        st = pychemia.structure_from_file(filename)
        self.assertEqual(st.nsites, 44)

    def test_from_file_2(self):
        """
        Test (pychemia.core.from_file)                              :
        """
        filename = 'tests/data/abinit_05/abinit.in'
        st = pychemia.structure_from_file(filename)
        self.assertEqual(st.nsites, 20)

    def test_from_file_3(self):
        """
        Test (pychemia.core.from_file)                              :
        """
        filename = 'tests/data/abinit_05/structure.json'
        st = pychemia.structure_from_file(filename)
        self.assertEqual(st.nsites, 20)
