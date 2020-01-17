import os
import shutil
import pychemia
import tempfile
import unittest


class MyTestCase(unittest.TestCase):

    def test_incar(self):
        """
        Test (pychemia.code.vasp) [INCAR parsing and writing]       :
        """
        print(os.getcwd())
        iv = pychemia.code.vasp.read_incar('tests/data/vasp_01/INCAR')
        self.assertEqual(len(iv), 12)
        self.assertEqual(iv.EDIFF, 1E-7)
        wf = tempfile.NamedTemporaryFile()
        iv.write(wf.name)
        wf.close()
        iv4dir = pychemia.code.vasp.read_incar('tests/data/vasp_01')
        self.assertEqual(iv, iv4dir)
        self.assertRaises(ValueError, pychemia.code.vasp.read_incar, 'tests/data')
        iv3 = pychemia.code.vasp.VaspInput(variables={'EDIFF': 1E-6})
        self.assertEqual(iv3['EDIFF'], 1E-6)
        iv = pychemia.code.vasp.read_incar('tests/data/vasp_02')
        iv.EDIFF *= 1.3
        td = tempfile.mkdtemp()
        pychemia.code.vasp.write_incar(iv, td)
        self.assertRaises(ValueError, iv.write_key, 'EDIF')
        shutil.rmtree(td)

    def test_bad_outcar(self):
        """
        Test (pychemia.code.vasp) [corrupted VASP OUTCAR]           :
        """
        vo = pychemia.code.vasp.VaspOutput('tests/data/vasp_04/OUTCAR')
        self.assertTrue(vo.is_finished)

    def test_encut_setup(self):
        """
        Test (pychemia.code.vasp) [ENCUT setup]                     :
        """
        iv = pychemia.code.vasp.read_incar('tests/data/vasp_06')
        iv.set_encut(ENCUT=1.2, POTCAR='tests/data/vasp_06/POTCAR')
        self.assertEqual(iv.ENCUT, 307)
        iv.set_rough_relaxation()
        self.assertEqual(iv.EDIFFG, -1E-2)
        iv.set_mit_settings()

    def test_vaspjob(self):
        """
        Test (pychemia.code.vasp) [VaspJob]                         :
        """
        td = tempfile.mkdtemp()
        st = pychemia.code.vasp.read_poscar('tests/data/vasp_06')
        kp = pychemia.code.vasp.read_kpoints('tests/data/vasp_06')
        self.assertEqual(kp.number_of_kpoints, 693)
        iv = pychemia.code.vasp.read_incar('tests/data/vasp_06')
        vj = pychemia.code.vasp.VaspJob(workdir=td,)
        vj.initialize(st, kpoints=kp)
        vj.set_input_variables(iv)
        vj.write_poscar()
        vj.write_kpoints()
        vj.write_incar()
        shutil.rmtree(td)

    def test_outcar(self):
        """
        Test (pychemia.code.vasp) [outcar]                          :
        """
        vo = pychemia.code.vasp.VaspOutput('tests/data/vasp_06/OUTCAR')
        self.assertEqual(vo.get_memory_used()['grid'], (1028.0, 'kBytes'))
        self.assertAlmostEqual(vo.to_dict['energy'], -19.67192646)
        print(vo)
        self.assertTrue(vo.has_forces_stress_energy())

    def test_poscar(self):
        """
        Test (pychemia.code.vasp) [poscar]                          :
        """
        # Temporal directory for outputs
        tmpdir = tempfile.mkdtemp()

        # Read a POSCAR by directory
        st = pychemia.code.vasp.read_poscar('tests/data/vasp_06')
        self.assertEqual(st.natom, 4)

        # Opening old format POSCAR without POTCAR
        with self.assertRaises(ValueError) as context:
            st = pychemia.code.vasp.read_poscar('tests/data/vasp_07/POSCAR')

        st = pychemia.code.vasp.read_poscar('tests/data/vasp_08/POSCAR_old')
        self.assertEqual(st.natom, 2)

        st = pychemia.code.vasp.read_poscar('tests/data/vasp_08/POSCAR_new')
        self.assertEqual(st.natom, 2)

        with self.assertRaises(ValueError) as context:
            pychemia.code.vasp.write_potcar(st, filepath=tmpdir + os.sep + 'POTCAR', basepsp='/no/existing/path')

        with self.assertRaises(ValueError) as context:
            pychemia.code.vasp.write_potcar(st, filepath=tmpdir + os.sep + 'POTCAR', basepsp='tests/data')

        cwd = os.getcwd()
        os.chdir('tests/data/vasp_07')
        st = pychemia.code.vasp.read_poscar('POSCAR_new')
        os.chdir(cwd)
        self.assertEqual(st.natom, 44)

        st = pychemia.code.vasp.read_poscar('tests/data/vasp_07/POSCAR_alt')

        pychemia.code.vasp.write_poscar(st, tmpdir + os.sep + 'POSCAR1')
        pychemia.code.vasp.write_poscar(st, tmpdir + os.sep + 'POSCAR2', direct=False)
        pychemia.code.vasp.write_poscar(st, tmpdir + os.sep + 'POSCAR3', newformat=False)

        st = pychemia.code.vasp.read_poscar(tmpdir + os.sep + 'POSCAR1')
        self.assertAlmostEqual(st.volume, 584.47161926043907)
        sym = pychemia.crystal.CrystalSymmetry(st)
        self.assertEqual(sym.symbol(), 'C2/c')
        st = pychemia.code.vasp.read_poscar(tmpdir + os.sep + 'POSCAR2')
        self.assertAlmostEqual(st.volume, 584.47161926043907)
        sym = pychemia.crystal.CrystalSymmetry(st)
        self.assertEqual(sym.symbol(), 'C2/c')
        st = pychemia.code.vasp.read_poscar(tmpdir + os.sep + 'POSCAR3')
        self.assertAlmostEqual(st.volume, 584.47161926043907)
        sym = pychemia.crystal.CrystalSymmetry(st)
        self.assertEqual(sym.symbol(), 'C2/c')
        pychemia.code.vasp.write_potcar(st, filepath=tmpdir + os.sep + 'POTCAR', basepsp='tests/data')
        pychemia.code.vasp.get_potcar_info(tmpdir + os.sep + 'POTCAR')

        shutil.rmtree(tmpdir)

