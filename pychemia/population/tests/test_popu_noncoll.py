#!/usr/bin/env python

import os
import pychemia
from pychemia.test import has_local_mongo


def test_popu_noncoll():
    """
    Tests (pychemia.population.NonCollinearMagMoms)              :
    """
    if not pychemia.HAS_PYMONGO:
        print('PyChemiaDB was disabled')
        return
    else:
        if not has_local_mongo():
            return

    source = 'pychemia/test/data/vasp_02'
    assert os.path.isfile(source + os.sep + 'INCAR')
    assert os.path.isfile(source + os.sep + 'POSCAR')
    popu = pychemia.population.NonCollinearMagMoms('test', source)
    popu.pcdb.clean()
    popu.random_population(16)

    assert len(popu) == 16
    popu.pcdb.clean()


if __name__ == '__main__':
    test_popu_noncoll()
