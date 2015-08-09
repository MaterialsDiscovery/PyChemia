#!/usr/bin/env python
__author__ = 'Guillermo Avendano-Franco'

import os
import pychemia


def test_popu_noncoll():
    """
    Testing PopulationNonColl           :
    """
    source = 'pychemia/test/data/vasp_02'
    assert os.path.isfile(source + os.sep + 'INCAR')
    assert os.path.isfile(source + os.sep + 'POSCAR')
    popu = pychemia.population.PopulationNonColl('test_PopulationNonColl', source)
    popu.random_population(32)

    assert len(popu) == 32
    popu.pcdb.clean()

if __name__ == '__main__':
    test_popu_noncoll()
