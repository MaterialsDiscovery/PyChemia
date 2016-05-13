import unittest
import numpy as np


class SearcherTest(unittest.TestCase):
    def test_harmony(self):
        """
        Tests (pychemia.searcher.harmony)                            :
        """
        import pychemia

        pychemia.pcm_log.debug('HarmonySearch')
        popu = pychemia.population.EuclideanPopulation(pychemia.utils.metaheuristics.Sphere.function, 3, [-1, 1])
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
