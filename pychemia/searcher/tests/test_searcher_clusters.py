import unittest

class SearcherTest(unittest.TestCase):

    def test_harmony(self):
        """
        Tests (pychemia.searcher.harmony) with LJ Clusters           :
        """
        import pychemia

        pychemia.pcm_log.debug('HarmonySearch')
        popu = pychemia.population.LJCluster('LJ', composition='Xe17')
        searcher = pychemia.searcher.HarmonySearch(popu, generation_size=16, stabilization_limit=5)
        searcher.run()
        popu.pcdb.clean()
        searcher = pychemia.searcher.HarmonySearch(popu, generation_size=16, stabilization_limit=5)
        searcher.run()

    def test_swarm(self):
        """
        Tests (pychemia.searcher.swarm) with LJ Clusters             :
        """
        import pychemia

        pychemia.pcm_log.debug('ParticleSwarm')
        popu = pychemia.population.LJCluster('LJ', composition='Xe17')
        searcher = pychemia.searcher.ParticleSwarm(popu, generation_size=16, stabilization_limit=5)
        searcher.run()
        popu.pcdb.clean()
        searcher = pychemia.searcher.ParticleSwarm(popu, generation_size=16, stabilization_limit=5)
        searcher.run()

    def test_firefly(self):
        """
        Tests (pychemia.searcher.firefly) with LJ Clusters           :
        """
        import pychemia

        pychemia.pcm_log.debug('FireFly')
        popu = pychemia.population.LJCluster('LJ', composition='Xe17')
        searcher = pychemia.searcher.FireFly(popu, generation_size=16, stabilization_limit=5)
        searcher.run()
        popu.pcdb.clean()
        searcher = pychemia.searcher.FireFly(popu,
                                             {'delta': 0.1,'gamma': 0.1, 'beta0': 0.8, 'alpha0': 0, 'multi_move': True},
                                             generation_size=16, stabilization_limit=5)
        searcher.run()

    def test_genetic(self):
        """

        Tests (pychemia.searcher.genetic) with LJ Clusters           :
        """
        import pychemia

        pychemia.pcm_log.debug('GeneticAlgorithm')
        popu = pychemia.population.LJCluster('LJ', composition='Xe17')
        searcher = pychemia.searcher.GeneticAlgorithm(popu, generation_size=16, stabilization_limit=5)
        searcher.run()
        popu.pcdb.clean()
        searcher = pychemia.searcher.GeneticAlgorithm(popu, generation_size=16, stabilization_limit=5)
        searcher.run()
