import pychemia
import unittest


class KPointsTest(unittest.TestCase):
    def test_optimized_grid(self):
        """
        Tests (pychemia.crystal.kpoints)                             :
        """
        st = pychemia.code.vasp.read_poscar('./data/SbBi/POSCAR')
        kp = pychemia.crystal.KPoints.optimized_grid(st.lattice)
        assert (kp.grid == [11, 11, 11])
        assert kp.get_density_of_kpoints(st.lattice) == 93126.888547945025
        kpdict = kp.to_dict
        kp2 = kp.from_dict(kpdict)
        assert kp.grid == kp2.grid
        kp = pychemia.crystal.KPoints(kmode='path',
                                      kvertices=[[0, 0, 0], [0.5, 0.5, 0], [0.5, 0.5, 0.5], [0, 0, 0]],
                                      intermediates=10)
        assert kp.number_of_kpoints == 40
        kp = pychemia.crystal.KPoints(kmode='reduced')
        kp.add_kpt([0, 0, 0], weight=0.5)
        assert kp.nkpt == 1
        kp.add_kpt([0.5, 0.5, 0.5], weight=0.5)
        print(kp)
        assert kp.number_of_kpoints == 2

        kp.set_kpoints_list([[0.0, 0.0, 0.0], [0.0, 0.5, 0.5], [0.5, 0.0, 0.5], [0.5, 0.5, 0.0], ])
        assert kp.number_of_kpoints == 4
