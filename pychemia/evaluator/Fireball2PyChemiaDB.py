from __future__ import print_function
import os
import socket
import time
import pychemia
from multiprocessing import Process
from pychemia import pcm_log, HAS_PYMONGO

if HAS_PYMONGO:
    from pychemia.db import get_database


class FireballCollector:

    def __init__(self, db_settings, dbname, source_dir='/fireball-archive/data-mining', nconcurrent=10, sleeping_time=2):
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


    def worker(self, db_settings, path):

        pcdb = get_database(db_settings)
        files = os.listdir(path)
        properties = {}
        geo = None

        if 'fireball.in' in files:
            score = 'input'
            try:
                invars = pychemia.code.fireball.read_fireball_in(path + os.sep + 'fireball.in')
                properties['input'] = invars
                properties['path'] = path
            except:
                print('Bad fireball.in')
                geo=None

            if 'param.dat' in files:
                try:
                    score += '+param'
                    param = pychemia.code.fireball.read_param(path + os.sep + 'param.dat')
                    properties['param_dat'] = param
                except:
                    print('Bad param.dat')
                    geo=none

            if 'eigen.dat' in files:
                try:
                    score += '+eigen'
                    eigen = pychemia.code.fireball.read_eigen(path + os.sep + 'eigen.dat')
                    properties['eigen_dat'] = eigen
                except:
                    print('bad eigen.dat')
                    geo=None

            if 'answer.bas' in files:
                try:
                    score += '+answer'
                    geo = pychemia.code.fireball.read_geometry_bas(path + os.sep + 'answer.bas')
                    periodic = False
                except:
                    print('Bad answer.bas')
                    geo=None

            for i in files:
                if i in self.outputs and os.path.isfile(path + os.sep + i):
                    rf = open(path + os.sep + i)
                    if "Welcome to FIREBALL" in rf.readline():
                        try:
                            output = pychemia.code.fireball.read_final_fireball_relax(path + os.sep + i)
                            properties['output']=output
                        except:
                            print('Bad %s' % i)
                    rf.close()
                if i[-3:]=='lvs':
                    try:
                        cell=pychemia.code.fireball.read_lvs(path + os.sep + i)
                        periodic=True
                    except:
                        print('Bad %s' % i)

        if geo is not None:
            print('DB: %d entries, Path: %s' % (pcdb.entries.count(), path))
            if periodic:
                st = pychemia.Structure(symbols=geo.symbols, positions=geo.positions, cell=cell)
            else:
                st = pychemia.Structure(symbols=geo.symbols, positions=geo.positions, periodicity=False)

            pcdb.insert(st, properties)
        return properties


    def process_directory(self, path, ret):

        files = os.listdir(path)
        score = False
        if os.path.exists(path+os.sep+'fireball.in') and os.path.exists(path+os.sep+'answer.bas'):
            for i in files:
                if i in self.outputs and os.path.isfile(path + os.sep + i):
                    rf = open(path + os.sep + i)
                    if "Welcome to FIREBALL" in rf.readline():
                        score = True
                    rf.close()

        if score:
            print(path)
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
        if not os.path.isfile('fireball_directories.txt'):
            self.process_directory(self.source_dir, ret)
            wf=open('fireball_directories.txt', 'w')
            for i in ret:
                wf.write("%s\n" % i)
            wf.close()
        else:
            rf=open('fireball_directories.txt')
            ret = [ x[:-1] for x in rf.readlines()]
            rf.close()

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
            
        to_evaluate = self.get_list_candidates()
        
        # Main loop looking permanently for candidates for evaluation
        while True:

            index = 0
            currently_evaluating = 0
            for j in range(self.nconcurrent):
                if procs[j] is not None and procs[j].is_alive():
                    currently_evaluating += 1
            print('Candidates to evaluate: %d  Candidates in evaluation: %d' % (len(to_evaluate), currently_evaluating))

            db_settings = dict(self.db_settings)
            # The first component of each pair in to_evaluate is the name of the database
            db_settings['name'] = self.dbname

            while index < len(to_evaluate):

                pcdb = get_database(db_settings)

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
                if not self.is_evaluated(pcdb, to_evaluate[index]):
                    ids_running[slot] = to_evaluate[index]
                    # This is the actual call to the worker, it must be a function with 4 arguments:
                    # The database settings, the entry identifier, the working directory and arguments for the worker
                    print('[Slot: %d] Evaluating: %s' % (slot,to_evaluate[index]))
                    procs[slot] = Process(target=self.worker, args=(db_settings, to_evaluate[index]))
                    procs[slot].start()
                    time.sleep(0.5)
                else:
                    print('[Slot: %d] Evaluated: %s' % (slot, to_evaluate[index]))
                index += 1
            time.sleep(self.sleeping_time)
