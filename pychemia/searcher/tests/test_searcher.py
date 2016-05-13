#!/usr/bin/env python

import unittest
import numpy as np
import time
import multiprocessing
# import logging

# logging.basicConfig(level=logging.DEBUG)
import pychemia
from pychemia import HAS_PYMONGO, HAS_GRIDFS


def evaluator():
    from pychemia.evaluator import cluster_launcher
    dbsettings = {'name': 'test'}
    cluster_launcher(dbsettings, 4)


def searcher():
    popu = pychemia.population.LJCluster('test', 'Ar13', refine=False, minimal_density=40)
    popu.pcdb.clean()
    hs = pychemia.searcher.HarmonySearch(popu, generation_size=8, stabilization_limit=3)
    hs.run()


def notest_searcher():
    """
    Testing HarmonySearch               :
    """
    if not HAS_PYMONGO or not HAS_GRIDFS:
        print('PyChemiaQueue was disabled')
        return

    p1 = multiprocessing.Process(target=searcher)
    p2 = multiprocessing.Process(target=evaluator)

    p1.start()
    time.sleep(10)
    p2.start()


class SearcherTest(unittest.TestCase):
    def notest_harmony(self):
        """
        Tests (pychemia.searcher.harmony)                            :
        """
        import pychemia

        pychemia.pcm_log.debug('HarmonySearch')
        popu = pychemia.population.EuclideanPopulation(pychemia.utils.metaheuristics.Sphere.function, 3, [-1, 1],
                                                       local_minimization=True)
        mini = pychemia.utils.metaheuristics.Sphere().minimum(3)
        hs = pychemia.searcher.HarmonySearch(popu)
        hs.run()
        assert np.linalg.norm(np.array(hs.population.db[hs.population.best_candidate]['x']) - mini) < 0.2

    def notest_swarm(self):
        """
        Tests (pychemia.searcher.swarm)                              :
        """
        import pychemia

        pychemia.pcm_log.debug('ParticleSwarm')
        popu = pychemia.population.EuclideanPopulation(pychemia.utils.metaheuristics.Sphere.function, 3, [-1, 1])
        mini = pychemia.utils.metaheuristics.Sphere().minimum(3)
        hs = pychemia.searcher.ParticleSwarm(popu)
        hs.run()
        assert np.linalg.norm(np.array(hs.population.db[hs.population.best_candidate]['x']) - mini) < 0.2

    def notest_firefly(self):
        """
        Tests (pychemia.searcher.firefly)                            :
        """
        import pychemia

        pychemia.pcm_log.debug('FireFly')
        popu = pychemia.population.EuclideanPopulation(pychemia.utils.metaheuristics.Sphere.function, 3, [-1, 1])
        mini = pychemia.utils.metaheuristics.Sphere().minimum(3)
        hs = pychemia.searcher.FireFly(popu)
        hs.run()
        assert np.linalg.norm(np.array(hs.population.db[hs.population.best_candidate]['x']) - mini) < 0.2

    def notest_genetic(self):
        """
        Tests (pychemia.searcher.genetic)                            :
        """
        import pychemia

        pychemia.pcm_log.debug('GeneticAlgorithm')
        popu = pychemia.population.EuclideanPopulation(pychemia.utils.metaheuristics.Sphere.function, 3, [-1, 1])
        mini = pychemia.utils.metaheuristics.Sphere().minimum(3)
        hs = pychemia.searcher.GeneticAlgorithm(popu)
        hs.run()
        assert np.linalg.norm(np.array(hs.population.db[hs.population.best_candidate]['x']) - mini) < 0.2


if __name__ == '__main__':
    test_searcher()
