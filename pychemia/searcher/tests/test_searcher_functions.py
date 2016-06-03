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

    def test_harmony(self):
        """
        Tests (pychemia.searcher.harmony)                            :
        """
        import pychemia

        pychemia.pcm_log.debug('HarmonySearch')
        mini = pychemia.utils.metaheuristics.Sphere().minimum(3)
        popu = pychemia.population.RealFunction(pychemia.utils.metaheuristics.Sphere.function, 3, [-1, 1],
                                                local_minimization=False)
        searcher = pychemia.searcher.HarmonySearch(popu, generation_size=16, stabilization_limit=5)
        searcher.run()
        popu = pychemia.population.RealFunction(pychemia.utils.metaheuristics.Sphere.function, 3, [-1, 1],
                                                local_minimization=True)

        searcher = pychemia.searcher.HarmonySearch(popu, generation_size=16, stabilization_limit=5)
        searcher.run()
        assert np.linalg.norm(np.array(searcher.population.db[searcher.population.best_candidate]['x']) - mini) < 0.2

    def test_swarm(self):
        """
        Tests (pychemia.searcher.swarm)                              :
        """
        import pychemia

        pychemia.pcm_log.debug('ParticleSwarm')
        mini = pychemia.utils.metaheuristics.Sphere().minimum(3)
        popu = pychemia.population.RealFunction(pychemia.utils.metaheuristics.Sphere.function, 3, [-1, 1],
                                                local_minimization=False)
        searcher = pychemia.searcher.ParticleSwarm(popu, generation_size=16, stabilization_limit=5)
        searcher.run()
        popu = pychemia.population.RealFunction(pychemia.utils.metaheuristics.Sphere.function, 3, [-1, 1],
                                                local_minimization=True)
        searcher = pychemia.searcher.ParticleSwarm(popu, generation_size=16, stabilization_limit=5)
        searcher.run()
        assert np.linalg.norm(np.array(searcher.population.db[searcher.population.best_candidate]['x']) - mini) < 0.2

    def test_firefly(self):
        """
        Tests (pychemia.searcher.firefly)                            :
        """
        import pychemia

        pychemia.pcm_log.debug('FireFly')
        mini = pychemia.utils.metaheuristics.Sphere().minimum(3)
        popu = pychemia.population.RealFunction(pychemia.utils.metaheuristics.Sphere.function, 3, [-1, 1],
                                                local_minimization=False)
        searcher = pychemia.searcher.FireFly(popu, generation_size=16, stabilization_limit=5)
        searcher.run()
        popu = pychemia.population.RealFunction(pychemia.utils.metaheuristics.Sphere.function, 3, [-1, 1],
                                                local_minimization=True)
        searcher = pychemia.searcher.FireFly(popu,
                                             {'delta': 0.1,'gamma': 0.1, 'beta0': 0.8, 'alpha0': 0, 'multi_move': True},
                                             generation_size=16, stabilization_limit=5)
        searcher.run()
        assert np.linalg.norm(np.array(searcher.population.db[searcher.population.best_candidate]['x']) - mini) < 0.2

    def test_genetic(self):
        """

        Tests (pychemia.searcher.genetic)                            :
        """
        import pychemia

        pychemia.pcm_log.debug('GeneticAlgorithm')
        mini = pychemia.utils.metaheuristics.Sphere().minimum(3)
        popu = pychemia.population.RealFunction(pychemia.utils.metaheuristics.Sphere.function, 3, [-1, 1],
                                                local_minimization=False)
        searcher = pychemia.searcher.GeneticAlgorithm(popu, generation_size=16, stabilization_limit=5)
        searcher.run()
        popu = pychemia.population.RealFunction(pychemia.utils.metaheuristics.Sphere.function, 3, [-1, 1],
                                                local_minimization=True)
        searcher = pychemia.searcher.GeneticAlgorithm(popu, generation_size=16, stabilization_limit=5)
        searcher.run()
        assert np.linalg.norm(np.array(searcher.population.db[searcher.population.best_candidate]['x']) - mini) < 0.2

