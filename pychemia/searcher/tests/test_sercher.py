#!/usr/bin/env python

import time
import multiprocessing
# import logging

#logging.basicConfig(level=logging.DEBUG)
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


if __name__ == '__main__':
    test_searcher()
