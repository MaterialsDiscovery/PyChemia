import unittest
import time
import multiprocessing
from pychemia.db import has_connection
from pychemia import pcm_log
from pychemia.searcher import HarmonySearch, FireFly, GeneticAlgorithm, ParticleSwarm
from pychemia.population import LJCluster


def evaluator():
    from pychemia.evaluator import cluster_launcher
    dbsettings = {'name': 'test'}
    cluster_launcher(dbsettings, 4)


def lj_searcher():
    popu = LJCluster('test', 'Ar13', refine=False, minimal_density=40)
    popu.pcdb.clean()
    hs = HarmonySearch(popu, generation_size=8, stabilization_limit=3)
    hs.run()


def notest_searcher():
    """
    Test (pychemia.searcher) with LJ Clusters                   :
    """
    if not has_connection():
        return

    p1 = multiprocessing.Process(target=lj_searcher)
    p2 = multiprocessing.Process(target=evaluator)

    p1.start()
    time.sleep(10)
    p2.start()


class SearcherTest(unittest.TestCase):

    def test_firefly(self):
        """
        Test (pychemia.searcher.firefly) with LJ Clusters           :
        """
        if not has_connection():
            return
        pcm_log.debug('FireFly')
        popu = LJCluster('test', composition='Xe13', refine=True, direct_evaluation=True)
        popu.pcdb.clean()
        searcher = FireFly(popu, generation_size=8, stabilization_limit=3)
        searcher.run()
        popu.pcdb.clean()
        searcher = FireFly(popu, {'delta': 0.1, 'gamma': 0.1, 'beta0': 0.8, 'alpha0': 0, 'multi_move': True},
                           generation_size=8, stabilization_limit=3)
        searcher.run()
        popu.pcdb.clean()

    def test_genetic(self):
        """

        Test (pychemia.searcher.genetic) with LJ Clusters           :
        """
        if not has_connection():
            return
        pcm_log.debug('GeneticAlgorithm')
        popu = LJCluster('test', composition='Xe13', refine=False, direct_evaluation=True)
        popu.pcdb.clean()
        searcher = GeneticAlgorithm(popu, generation_size=8, stabilization_limit=3)
        searcher.run()
        popu.pcdb.clean()

    def test_harmony(self):
        """
        Test (pychemia.searcher.harmony) with LJ Clusters           :
        """
        if not has_connection():
            return
        pcm_log.debug('HarmonySearch')
        popu = LJCluster('test', composition='Xe13', refine=False, direct_evaluation=True)
        popu.pcdb.clean()
        searcher = HarmonySearch(popu, generation_size=8, stabilization_limit=3)
        searcher.run()
        popu.pcdb.clean()

    def test_swarm(self):
        """
        Test (pychemia.searcher.swarm) with LJ Clusters             :
        """
        if not has_connection():
            return
        pcm_log.debug('ParticleSwarm')
        popu = LJCluster('test', composition='Xe13', refine=False, direct_evaluation=True)
        popu.pcdb.clean()
        searcher = ParticleSwarm(popu, generation_size=8, stabilization_limit=3)
        searcher.run()
        popu.pcdb.clean()

