#!/usr/bin/env python

import os
import pychemia
from pychemia import HAS_PYMONGO

def test_popu_noncoll():
    """
    Tests (pychemia.population.PopulationNonColl)                :
    """
    if not HAS_PYMONGO:
        print('PyChemiaDB was disabled')
        return

    try:
        import pymongo
        maxSevSelDelay = 1
        client = pymongo.MongoClient("localhost", serverSelectionTimeoutMS=maxSevSelDelay)
        client.server_info()  # force connection on a request as the
        # connect=True parameter of MongoClient seems
        # to be useless here
    except pymongo.errors.ServerSelectionTimeoutError as err:
        # do whatever you need
        print(err)
        return

    source = 'pychemia/test/data/vasp_02'
    assert os.path.isfile(source + os.sep + 'INCAR')
    assert os.path.isfile(source + os.sep + 'POSCAR')
    popu = pychemia.population.PopulationNonColl('test_PopulationNonColl', source)
    popu.random_population(32)

    assert len(popu) == 32
    popu.pcdb.clean()


if __name__ == '__main__':
    test_popu_noncoll()
