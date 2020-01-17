import pychemia
import unittest
import numpy as np


class CrystalSymmetryTest(unittest.TestCase):
    def test_optimized_grid(self):
        """
        Test (pychemia.crystal.symmetry)                            :
        """
        from pychemia import pcm_log
        pcm_log.debug("CrystalSymmetryTest")
        st = pychemia.code.vasp.read_poscar('tests/data/SbBi/POSCAR')
        cs = pychemia.crystal.CrystalSymmetry(st)
        assert cs.number() == 160
        assert cs.symbol() == u'R3m'
        pr = cs.find_primitive()
        assert np.abs(st.volume - pr.volume) < 1E-6
        assert cs.crystal_system() == u'Trigonal'
        ss = cs.symmetrize()
        assert np.abs(st.volume - ss.volume) < 1E-6
