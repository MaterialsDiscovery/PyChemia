from __future__ import unicode_literals
import unittest
# import unittest2 as unittest
import os
from pychemia import HAS_PYMONGO

if HAS_PYMONGO:
    from pychemia.population import LJCluster, StructurePopulation, PopulationDFTU, PopulationNonColl, \
        EuclideanPopulation


def has_connection():
    if not HAS_PYMONGO:
        return False
    import pymongo
    try:
        maxSevSelDelay = 1
        client = pymongo.MongoClient("localhost", serverSelectionTimeoutMS=maxSevSelDelay)
        client.server_info()  # force connection on a request as the
        # connect=True parameter of MongoClient seems
        # to be useless here
        return True
    except pymongo.errors.ServerSelectionTimeoutError as err:
        # do whatever you need
        print(err)
        return False


def funx2(x):
    return x ** 2


class PopulationTest(unittest.TestCase):
    def test_ljcluster(self):
        """
        Tests (pychemia.population.LJCluster)                        :
        """
        if not has_connection():
            return
        popu = LJCluster('test', 'Ne4')
        popu.add_random()
        popu.add_random()
        popu.pcdb.clean()

    def test_structure(self):
        """
        Tests (pychemia.population.StructurePopulation)              :
        """
        if not has_connection():
            return
        popu = StructurePopulation('test', 'NaCl')
        popu.add_random()
        popu.add_random()
        popu.pcdb.clean()

    def test_noncoll(self):
        """
        Tests (pychemia.population.PopulationNonColl)                :
        """
        if not has_connection():
            return
        popu = PopulationNonColl('test', source_dir='pychemia/test/data/vasp_02')
        popu.add_random()
        popu.add_random()
        popu.pcdb.clean()

    def notest_dftu(self):
        """
        Tests (pychemia.population.PopulationDFTU)                   :
        """
        if not has_connection():
            return
        print(os.getcwd())
        popu = PopulationDFTU('test', './data/abinit_01/abinit.in', oxidations=[-1, 1])
        popu.add_random()
        popu.add_random()
        popu.pcdb.clean()

    def test_euclidean(self):
        """
        Tests (pychemia.population.EuclideanPopulation)              :
        """
        if not has_connection():
            return
        popu = EuclideanPopulation(funx2, 2, [-1, 1])
        popu.add_random()
        popu.add_random()


if __name__ == "__main__":
    unittest.main()
