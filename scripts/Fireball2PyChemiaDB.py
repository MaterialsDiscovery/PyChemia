
import os
import socket
import time
import pychemia
from multiprocessing import Process
from pychemia import pcm_log, HAS_PYMONGO

if HAS_PYMONGO:
    from pychemia.db import get_database


class FireballCollector:

    def __init__(self, db_settings, dbname, source_dir='/fireball-archive/data-mining', nconcurrent=1,
                 sleeping_time=120):
        """
        FireballCollector is a class to manage the execution of a function 'worker' for entries on a list of PyChemiaDB
         databases.
        The execution of each worker occurs directly on the machine executing the code, and the class controls the
        number of concurrent execution using the variable 'nconcurrent'

        :param db_settings: Common database settings for all the databases that will be processed by this evaluator.
                            All databases should share common settings like server name, user and password, ssl and
                            replicaset settings.
        :param dbname: List of database names, the common settings are stored in db_settings
        :param nconcurrent: Number of concurrent executions allowed. Their number usually should be the number of cores
                            available on the machine. (default: 1)
        :param sleeping_time: Time in seconds before try to search for new candidates to evaluate (default: 120 seconds)
        """
        self.db_settings = db_settings
        self.dbname = dbname
        self.source_dir = source_dir
        self.nconcurrent = nconcurrent
        self.sleeping_time = sleeping_time

        rf = open(self.source_dir+os.sep+'outputs.txt')
        self.outputs = [x.strip() for x in rf.readlines()]

    def worker(self, pcdb, path):

        files = os.listdir(path)
        properties = {}
        geo = None

        if 'fireball.in' in files:
            score = 'input'
            invars = pychemia.code.fireball.read_fireball_in(path + os.sep + 'fireball.in')
            properties['input'] = invars
            if 'param.dat' in files:
                score += '+param'
                param = pychemia.code.fireball.read_param(path + os.sep + 'param.dat')
                properties['param_dat'] = param
            if 'eigen.dat' in files:
                score += '+eigen'
                eigen = pychemia.code.fireball.read_eigen(path + os.sep + 'eigen.dat')
                properties['eigen_dat'] = eigen
            if 'answer.bas' in files:
                score += '+answer'
                geo = pychemia.code.fireball.read_geometry_bas(path + os.sep + 'answer.bas')

            for i in files:
                if i in self.outputs and os.path.isfile(path + os.sep + i):
                    rf = open(path + os.sep + i)
                    if "Welcome to FIREBALL" in rf.readline():
                        outputs = pychemia.code.fireball.read_fireball_stdout(path + os.sep + i)
                        properties['output'] = outputs
                    rf.close()

        if geo is not None:
            pcdb.insert(geo, properties)

    def process_directory(self, path, ret):

        files = os.listdir(path)
        score = True
        if 'fireball.in' in files and 'answer.bas' in files:
            for i in files:
                if os.path.isfile(path + os.sep + i):
                    rf = open(path + os.sep + i)
                    if "Welcome to FIREBALL" in rf.readline():
                        score = True
                    rf.close()

        if score:
            ret.append(path)

        for i in files:
            if os.path.isdir(path + os.sep + i):
                self.process_directory(path + os.sep + i, ret)

    def get_list_candidates(self):
        """
        Scan all databases looking for candidates for evaluation

        :return: A list of pairs, each pair contains the name of the database
                    and candidate identifier.
        """

        ret = []
        self.process_directory(self.source_dir, ret)
        return ret

    def is_evaluated(self, pcdb, path):
        if pcdb.entries.find({'properties.path': path}).count() > 0:
            return True
        else:
            return False

    def run(self):
        """
        Continuously search for suitable candidates to evaluation among a list of databases.

        :return:
        """

        procs = []
        ids_running = []

        # Creating a list to store the 'nconcurrent' running jobs
        for i in range(self.nconcurrent):
            procs.append(None)
            ids_running.append(None)

        # Main loop looking permanently for candidates for evaluation
        while True:

            to_evaluate = self.get_list_candidates()

            index = 0
            currently_evaluating = 0
            for j in range(self.nconcurrent):
                if procs[j] is not None and procs[j].is_alive():
                    currently_evaluating += 1
            print('Candidates to evaluate: %d  Candidates in evaluation: %d' % (len(to_evaluate), currently_evaluating))

            while index < len(to_evaluate):

                db_settings = dict(self.db_settings)
                # The first component of each pair in to_evaluate is the name of the database
                db_settings['name'] = self.dbname
                pcdb = get_database(db_settings)

                slot = None
                while True:
                    for j in range(self.nconcurrent):
                        if procs[j] is None or not procs[j].is_alive():
                            slot = j
                    if slot is None:
                        time.sleep(self.sleeping_time)
                    else:
                        break

                # The function is_evaluated needs two arguments, the database object and entry identifier and
                # must return a boolean to decide if the candidate should be evaluated.
                if self.is_evaluated(pcdb, to_evaluate[index]):
                    ids_running[slot] = to_evaluate[index]

                    # This is the actual call to the worker, it must be a function with 4 arguments:
                    # The database settings, the entry identifier, the working directory and arguments for the worker
                    procs[j] = Process(target=self.worker, args=(db_settings, to_evaluate[index]))
                    procs[j].start()
                else:
                    pcm_log.debug('Not evaluable: %s' % to_evaluate[index])
                index += 1
            time.sleep(self.sleeping_time)
