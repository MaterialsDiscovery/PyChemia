
import os
import socket
import time
from multiprocessing import Process
from pychemia import pcm_log, HAS_PYMONGO

if HAS_PYMONGO:
    from pychemia.db import get_database


class DirectEvaluator:

    def __init__(self, db_settings, dbnames, source_dir, is_evaluated, worker, worker_args=None, nconcurrent=1,
                 evaluate_failed=False, evaluate_all=False, sleeping_time=120):
        """
        DirectEvaluator is a class to manage the execution of a function 'worker' for entries on a list of PyChemiaDB
         databases.
        The execution of each worker occurs directly on the machine executing the code, and the class controls the
        number of concurrent execution using the variable 'nconcurrent'

        :param db_settings: Common database settings for all the databases that will be processed by this evaluator.
                            All databases should share common settings like server name, user and password, ssl and
                            replicaset settings.
        :param dbnames: List of database names, the common settings are stored in db_settings
        :param source_dir: Path to directory where the candidates will be evaluated.
        :param is_evaluated: Function to decide if one entry is evaluable
        :param worker: Function applied to each candidate on the databases that are not considered evaluated.
        :param worker_args: Arguments for the worker function. (Usually a dictionary that depends on the arguments
                            needed by worker)
        :param nconcurrent: Number of concurrent executions allowed. Their number usually should be the number of cores
                            available on the machine. (default: 1)
        :param evaluate_failed: Boolean to decide if candidates marked as failed should be evaluated again.
        :param evaluate_all: Boolean to decide is all candidates are evaluated regardless of the outcome of the function
                is_evaluated. (default: False)
        :param sleeping_time: Time in seconds before try to search for new candidates to evaluate (default: 120 seconds)
        """
        self.db_settings = db_settings
        self.dbnames = dbnames
        self.source_dir = source_dir
        self.is_evaluated = is_evaluated
        self.worker = worker
        self.worker_args = worker_args
        self.nconcurrent = nconcurrent
        self.evaluate_failed = evaluate_failed
        self.evaluate_all = evaluate_all
        self.sleeping_time = sleeping_time

    def unlock_all(self):
        """
        Checking all databases and unlocking all entries.

        :return: None
        """
        for idb in self.dbnames:
            db_settings = dict(self.db_settings)
            db_settings['name'] = idb
            pcdb = get_database(db_settings)
            print('Database contains: %d entries' % pcdb.entries.count())

            for entry in pcdb.entries.find({}):
                pcdb.unlock(entry['_id'], name=socket.gethostname())
            print('Number of entries: ', pcdb.entries.count())

    def get_list_candidates(self):
        """
        Scan all databases looking for candidates for evaluation

        :return: A list of pairs, each pair contains the name of the database
                    and candidate identifier.
        """
        ret = []
        for idb in self.dbnames:
            print(idb)
            db_settings = dict(self.db_settings)
            db_settings['name'] = idb
            pcdb = get_database(db_settings)

            for entry in pcdb.entries.find({}, {'_id': 1}):
                entry_id = entry['_id']
                if not self.is_evaluated(pcdb, entry_id, self.worker_args) or self.evaluate_all:
                    pcm_log.debug('Adding entry %s from db %s' % (str(entry_id), pcdb.name))
                    ret.append([idb, entry_id])
        print('Found %d entries to evaluate' % len(ret))
        return ret

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

        self.unlock_all()

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
                dbname = to_evaluate[index][0]
                db_settings['name'] = dbname
                pcdb = get_database(db_settings)
                # The second component of each pair in to_evaluate is the entry_id
                entry_id = to_evaluate[index][1]

                for j in range(self.nconcurrent):
                    if procs[j] is None or not procs[j].is_alive():
                        ids_running[j] = None

                if entry_id in ids_running:
                    print('Already executing: %s' % entry_id)
                    index += 1
                    continue
                else:
                    print('DB: %10s Entry: %s' % (dbname, entry_id))

                if not os.path.exists(self.source_dir + os.sep + dbname):
                    os.mkdir(self.source_dir + os.sep + dbname)

                slot = None
                while True:
                    for j in range(self.nconcurrent):
                        if procs[j] is None or not procs[j].is_alive():
                            slot = j
                            break
                    if slot is None:
                        time.sleep(self.sleeping_time)
                    else:
                        break

                # The function is_evaluated needs two arguments, the database object and entry identifier and
                # must return a boolean to decide if the candidate should be evaluated.
                if not self.is_evaluated(pcdb, entry_id, self.worker_args) or self.evaluate_all:
                    pcm_log.debug('Evaluable: %s:%s. Relaxing entry %d of %d Slot: %d' % (dbname,
                                                                                          str(entry_id),
                                                                                          index,
                                                                                          len(to_evaluate), slot))
                    ids_running[slot] = entry_id
                    workdir = self.source_dir + os.sep + dbname + os.sep + str(entry_id)
                    if not os.path.exists(workdir):
                        os.mkdir(workdir)
                    pcm_log.debug('Launching for %s id: %s' % (pcdb.name, str(entry_id)))

                    # This is the actual call to the worker, it must be a function with 4 arguments:
                    # The database settings, the entry identifier, the working directory and arguments for the worker
                    procs[slot] = Process(target=self.worker, args=(db_settings, entry_id, workdir, self.worker_args))
                    procs[slot].start()
                    time.sleep(1)
                else:
                    pcm_log.debug('Not evaluable: %s' % str(entry_id))
                index += 1
            time.sleep(self.sleeping_time)
