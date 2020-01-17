import unittest
import pychemia
from pychemia.utils.serializer import generic_serializer
import numpy as np
import tempfile
import shutil


class PopulationTest(unittest.TestCase):
    def test_orbital(self):
        """
        Test (pychemia.population.OrbitalDFTU)                      :
        """
        if not pychemia.db.has_connection():
            return

        pychemia_path = pychemia.__path__[0]
        abiinput = pychemia.code.abinit.AbinitInput('tests/data/abinit_dmatpawu/abinit.in')
        dmatpawu = np.array(abiinput['dmatpawu']).reshape(-1, 5, 5)
        params = pychemia.population.orbitaldftu.dmatpawu2params(dmatpawu, 5)
        dmatpawu_new = pychemia.population.orbitaldftu.params2dmatpawu(params)

        self.assertLess(np.min(dmatpawu-dmatpawu_new), 0.01)

        with self.assertRaises(ValueError) as context:
            pychemia.population.orbitaldftu.OrbitalDFTU('test', '/tmp/no_abinit.in')

        with self.assertRaises(ValueError) as context:
            pychemia.population.orbitaldftu.OrbitalDFTU('test', input_path='tests/data/abinit_01/abinit.in')

        popu = pychemia.population.orbitaldftu.OrbitalDFTU('test', input_path='tests/data/abinit_dmatpawu/abinit.in')
        popu.pcdb.clean()

        ea = [[-0.1, -0.2, -0.3, -0.4, -0.5, -0.6, -0.7, -0.8, -0.9, -1.0],
              [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0],
              [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
              [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]]
        params['euler_angles'] = ea
        params = generic_serializer(params)
        entry_id = popu.new_entry(params)

        ea = [[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0],
              [-0.1, -0.2, -0.3, -0.4, -0.5, -0.6, -0.7, -0.8, -0.9, -1.0],
              [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
              [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]]
        params['euler_angles'] = ea
        params = generic_serializer(params)
        entry_jd = popu.new_entry(params)

        popu.add_random()
        popu.random_population(16)
        print(popu)

        self.assertFalse(popu.is_evaluated(entry_id))

        for entry_id in popu.members:
            params = popu.get_correlation_params(entry_id, final=False)
            popu.set_final_results(entry_id, params, 0.0, 1E-13)

        self.assertTrue(popu.is_evaluated(entry_id))

        popu.get_duplicates(popu.members, tolerance=0.1)

        popu.cross([entry_id, entry_jd])

        entry_idm = popu.move_random(entry_id)
        popu.get_entry(entry_idm, {'properties': 1})

        entry_imj = popu.move(entry_id, entry_jd)
        popu.get_entry(entry_imj, {'properties': 1})

        # pd = popu.to_dict
        # popu.from_dict(pd)

        tmpdir = tempfile.mkdtemp()

        for i in popu.members:
            popu.prepare_folder(i, workdir=tmpdir, source_dir=pychemia_path + 'tests/data/abinit_dmatpawu')

        popu.pcdb.clean()
        shutil.rmtree(tmpdir)
