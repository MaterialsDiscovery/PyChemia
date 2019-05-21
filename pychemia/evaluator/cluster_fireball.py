
import os
import time
from multiprocessing import Pool, Process
import pychemia
from pychemia.utils.serializer import generic_serializer

__author__ = 'Guillermo Avendano-Franco'


def cluster_fb_worker(db_settings):
    while True:
        pcdb = pychemia.db.get_database(db_settings)
        population = pychemia.population.LJCluster(pcdb)

        entry = population.pcdb.db.pychemia_entries.find_one({'status.' + population.tag: True,
                                                              'status.lock': {'$exists': False},
                                                              'properties': {}}, {'_id': 1})
        if entry is not None:
            population.pcdb.lock(entry['_id'])
            structure = population.pcdb.get_structure(entry['_id'])

            fb = pychemia.code.fireball.FireBall(fdata_path='../Fdata')
            fb.initialize(structure, workdir=str(entry['_id']))
            fb.cluster_relaxation()
            fb.set_inputs()
            sp = fb.run()
            sp.wait()
            so = pychemia.code.fireball.read_fireball_stdout(fb.workdir + os.sep + 'fireball.log')
            forces = generic_serializer(so['forces'][-1])
            energy = so['energy'][-1]['ETOT']
            properties = {'forces': forces, 'energy': energy}
            structure = pychemia.code.fireball.read_geometry_bas(fb.workdir + os.sep + 'answer.bas')
            population.pcdb.update(entry['_id'], structure=structure, properties=properties)
            population.pcdb.unlock(entry['_id'])
        else:
            break


def cluster_fb_evaluator(db_settings, nparal):
    pcdb = pychemia.db.get_database(db_settings)
    population = pychemia.population.LJCluster(pcdb)
    print('Staring evaluator for ', population.name)
    while True:
        entry = population.pcdb.db.pychemia_entries.find_one({'status.' + population.tag: True,
                                                              'status.lock': {'$exists': False},
                                                              'properties': {}}, {'_id': 1})

        if entry is None:
            time.sleep(2)
            create_pool = False
        else:
            create_pool = True

        if create_pool:
            pool = Pool(processes=nparal)
            pool.map(cluster_fb_worker, nparal * [db_settings])
            pool.close()
            pool.join()


def cluster_fb_launcher(db_settings, nparal):
    p = Process(target=cluster_fb_evaluator, args=(db_settings, nparal))
    p.start()
    return p
