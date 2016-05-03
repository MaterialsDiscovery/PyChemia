#!/usr/bin/env python

import os
import pychemia
from pychemia import HAS_PYMONGO


def test_popu_noncoll():
    """
    Testing PopulationNonColl           :
    """
    if not HAS_PYMONGO:
        print('PyChemiaDB was disabled')
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
