import os
import pychemia
import tempfile
import unittest


class MyTestCase(unittest.TestCase):

    def test_incar(self):
        """
        Tests (pychemia.code.vasp) [INCAR parsing and writing]       :
        """
        print(os.getcwd())
        iv = pychemia.code.vasp.read_incar('pychemia/test/data/vasp_01/INCAR')
        self.assertEqual(len(iv), 12)
        self.assertEqual(iv.EDIFF, 1E-7)
        wf = tempfile.NamedTemporaryFile()
        iv.write(wf.name)
        wf.close()
        iv4dir = pychemia.code.vasp.read_incar('pychemia/test/data/vasp_01')
        self.assertEqual(iv, iv4dir)
        self.assertRaises(ValueError, pychemia.code.vasp.read_incar, 'pychemia/test/data')
        iv3 = pychemia.code.vasp.InputVariables(variables={'EDIFF': 1E-6})
        self.assertEqual(iv3['EDIFF'], 1E-6)
        iv = pychemia.code.vasp.read_incar('pychemia/test/data/vasp_02')
        iv.EDIFF *= 1.3
        td = tempfile.TemporaryDirectory()
        pychemia.code.vasp.write_incar(iv, td.name)
        self.assertRaises(ValueError, iv.write_key, 'EDIF')

    def test_bad_outcar(self):
        """
        Tests (pychemia.code.vasp) [corrupted VASP OUTCAR]           :
        """
        vo = pychemia.code.vasp.VaspOutput('pychemia/test/data/vasp_04/OUTCAR')
        self.assertTrue(vo.is_finished)

    def test_encut_setup(self):
        """
        Tests (pychemia.code.vasp) [ENCUT setup]                     :
        """
        iv = pychemia.code.vasp.read_incar('pychemia/test/data/vasp_06')
        iv.set_encut(ENCUT=1.2, POTCAR='pychemia/test/data/vasp_06/POTCAR')
        self.assertEqual(iv.ENCUT, 307)
        iv.set_rough_relaxation()
        self.assertEqual(iv.EDIFFG, -1E-2)
        iv.set_mit_settings()

    def test_vaspjob(self):
        """
        Tests (pychemia.code.vasp) [VaspJob]                         :
        """
        td = tempfile.TemporaryDirectory()
        st=pychemia.code.vasp.read_poscar('pychemia/test/data/vasp_06')
        kp=pychemia.code.vasp.read_kpoints('pychemia/test/data/vasp_06')
        self.assertEqual(kp.number_of_kpoints, 693)
        iv = pychemia.code.vasp.read_incar('pychemia/test/data/vasp_06')
        vj=pychemia.code.vasp.VaspJob()
        vj.initialize(st, workdir=td.name, kpoints=kp)
        vj.set_input_variables(iv)
        vj.write_poscar()
        vj.write_kpoints()
        vj.write_incar()

    def test_poscar(self):
        """
        Tests (pychemia.code.vasp) [VaspJob]                         :
        """
        st=pychemia.code.vasp.read_poscar('pychemia/test/data/vasp_06')
