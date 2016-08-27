from __future__ import unicode_literals
import unittest
# import unittest2 as unittest
import os
from pychemia import HAS_PYMONGO

if HAS_PYMONGO:
    from pychemia.population import LJCluster, RelaxStructures, OrbitalDFTU, NonCollinearMagMoms, \
        RealFunction


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

    def notest_structure(self):
        """
        Tests (pychemia.population.RelaxStructures)                  :
        """
        if not has_connection():
            return
        popu = RelaxStructures('test', 'NaCl')
        popu.add_random()
        popu.add_random()
        popu.pcdb.clean()

    def test_noncoll(self):
        """
        Tests (pychemia.population.NonCollinearMagMoms)              :
        """
        if not has_connection():
            return
        popu = NonCollinearMagMoms('test', source_dir='pychemia/test/data/vasp_02')
        popu.add_random()
        popu.add_random()
        popu.pcdb.clean()

    def notest_dftu(self):
        """
        Tests (pychemia.population.OrbitalDFTU)                      :
        """
        if not has_connection():
            return
        print(os.getcwd())
        popu = OrbitalDFTU('test', './data/abinit_05/abinit.in')
        popu.add_random()
        popu.add_random()
        popu.pcdb.clean()

    def test_euclidean(self):
        """
        Tests (pychemia.population.RealFunction)                     :
        """
        if not has_connection():
            return
        popu = RealFunction(funx2, 2, [-1, 1])
        popu.add_random()
        popu.add_random()


if __name__ == "__main__":
    unittest.main()
