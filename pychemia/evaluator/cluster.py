import time
from multiprocessing import Pool, Process
import pychemia

__author__ = 'Guillermo Avendano-Franco'


def cluster_worker(db_settings):
    while True:
        pcdb = pychemia.db.get_database(db_settings)
        population = pychemia.population.LJCluster(pcdb)

        entry = population.pcdb.db.pychemia_entries.find_one({'status.' + population.tag: True,
                                                              'status.lock': {'$exists': False},
                                                              'properties': {}}, {'_id': 1})
        if entry is not None:
            population.pcdb.lock(entry['_id'])
            structure, properties, energy = population.evaluate(entry['_id'])
            population.pcdb.update(entry['_id'], structure=structure, properties=properties)
            population.pcdb.unlock(entry['_id'])
        else:
            break


def cluster_evaluator(db_settings, nparal):
    pcdb = pychemia.db.get_database(db_settings)
    population = pychemia.population.LJCluster(pcdb)
    population.recover()
    print('Staring evaluator for %s' % population.name)
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
            pool.map(cluster_worker, nparal * [db_settings])
            pool.close()
            pool.join()


def cluster_launcher(db_settings, nparal):
    p = Process(target=cluster_evaluator, args=(db_settings, nparal))
    p.start()
    return p
