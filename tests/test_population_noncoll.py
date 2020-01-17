#!/usr/bin/env python

import os
import pychemia
from .local_mongo import has_local_mongo


def notest_popu_noncoll():
    """
    Test (pychemia.population.NonCollinearMagMoms)              :
    """
    if not pychemia.HAS_PYMONGO:
        print('PyChemiaDB was disabled')
        return
    else:
        if not has_local_mongo():
            return

    source = 'tests/data/vasp_02'
    assert os.path.isfile(source + os.sep + 'INCAR')
    assert os.path.isfile(source + os.sep + 'POSCAR')
    popu = pychemia.population.NonCollinearMagMoms('test', source)
    popu.pcdb.clean()
    popu.random_population(16)

    assert len(popu) == 16
    popu.pcdb.clean()


if __name__ == '__main__':
    test_popu_noncoll()
