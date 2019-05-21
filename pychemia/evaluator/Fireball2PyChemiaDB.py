
import os
import socket
import time
import pychemia
from multiprocessing import Process
from pychemia import pcm_log, HAS_PYMONGO

if HAS_PYMONGO:
    from pychemia.db import get_database


class FireballCollector:

    def __init__(self, dbname, source_dir=None, source_file=None, output_names=None, output_file=None,
                 db_settings=None, nconcurrent=1):
        """
        FireballCollector is a class collect Fireball executions into one PyChemia Database.

        The name of the database will be 'dbname' and executions will be processed from a give 'source_dir', a path
        to search for Fireball calculations or 'source_file', one file with one directory per line where Fireball
        executions are stored.
        The variable 'output_names' sets possible output files for Fireball. Typically Fireball output is the standard
        output and the user has to redirect the output to a file with an arbitrary name.

        The processing of directories happens in parallel using 'nconcurrent' python processes.
        By default the database is created on the local machine assuming no authorization or authentication, otherwise
        use the dictionary 'db_settings' to get more control on the database location and access.

        :param dbname:      (str) Name of the database that will contain the Fireball executions
        :param source_dir:  (str) Path for a directory that will be explored recursevely for suitable Fireball
                                  executions
        :param source_file: (str) Path to a file describing directories with Fireball executions, one per line
        :param output_names: (list) List of possible output filenames
        :param output_file: (str) Path to a file with possible output names, one per line
        :param db_settings: (dict) Extra parameters for locating and accessing the mongo server.
        :param nconcurrent: Number of concurrent python executions to fill the database (default: 1)
        """
        self.db_settings = db_settings
        self.dbname = dbname
        self.source_dir = source_dir
        self.source_file = source_file
        self.output_names = output_names
        self.output_file = output_file
        self.nconcurrent = nconcurrent
        self.sleeping_time = 2
        self.database = None

        if self.output_file is not None and os.path.exists(self.output_file):
            rf = open(self.output_file)
            self.output_names = [x.strip() for x in rf.readlines()]

        if self.source_dir is None and self.source_file is None:
            raise ValueError("One of the variables source_dir or source_file need to be defined")

        if self.db_settings is None:
            self.db_settings = {'name': self.dbname}
        else:
            self.db_settings['name'] = self.dbname

        self.db_settings = dict(self.db_settings)
        # The first component of each pair in to_evaluate is the name of the database
        self.db_settings['name'] = self.dbname

        self.database = get_database(self.db_settings)

    def worker(self, path):
        """
        From a given 'path', the worker will collect information from several files
        such as:

        fireball.in
        eigen.dat
        param.dat
        answer.bas
        output collected from the standard output of fireball.x


        After collecting data into dictionaries, the data is stored in a PyChemia database
        """

        pcdb = get_database(self.db_settings)
        files = os.listdir(path)
        properties = {}
        geo = None

        if os.path.lexists(path+os.sep+'fireball.in'):
            properties['path'] = path
            score = 'input'
            try:
                invars = pychemia.code.fireball.read_fireball_in(path + os.sep + 'fireball.in')
                properties['input'] = invars
            except:
                print('Bad fireball.in on %s' % path)

            if os.path.exists(path+os.sep+'param.dat'):
                try:
                    score += '+param'
                    param = pychemia.code.fireball.read_param(path + os.sep + 'param.dat')
                    properties['param_dat'] = param
                except:
                    print('Bad param.dat')

            if os.path.lexists(path+os.sep+'eigen.dat'):
                try:
                    score += '+eigen'
                    eigen = pychemia.code.fireball.read_eigen(path + os.sep + 'eigen.dat')
                    properties['eigen_dat'] = eigen
                except:
                    print('bad eigen.dat')

            if os.path.lexists(path+os.sep+'answer.bas'):
                try:
                    score += '+answer'
                    geo = pychemia.code.fireball.read_geometry_bas(path + os.sep + 'answer.bas')
                    periodic = False
                except:
                    print('Bad answer.bas')

            for ioutput in self.output_names:
                if os.path.isfile(path + os.sep + ioutput):
                    rf = open(path + os.sep + ioutput)
                    line = rf.readline()
                    rf.close()
                    if "Welcome to FIREBALL" in line:
                        try:
                            output = pychemia.code.fireball.read_final_fireball_relax(path + os.sep + ioutput)
                            properties['output'] = output
                            break
                        except:
                            print('Bad output %s on %s' % (ioutput, path))

            for ifile in files:
                if ifile[-3:] == 'lvs':
                    try:
                        cell = pychemia.code.fireball.read_lvs(path + os.sep + ifile)
                        periodic = True
                        lat = pychemia.crystral.Lattice(cell)
                        if lat.volume > 1E7:
                            print('Lattice too big, assuming non-periodic structure')
                            periodic = False
                    except:
                        print('Bad %s' % ifile)

        if geo is not None:
            print('DB: %d entries, Path: %s' % (pcdb.entries.count(), path))
            if periodic:
                st = pychemia.Structure(symbols=geo.symbols, positions=geo.positions, cell=cell)
            else:
                st = pychemia.Structure(symbols=geo.symbols, positions=geo.positions, periodicity=False)

            pcdb.insert(st, properties)
        return properties

    def process_directory(self, path, fireball_dirs):
        """
        Search for directories with fireball.in and answer.bas and add them to the
        current list of 'fireball_dirs'

        The search is recursive to subdirectories.
        """
        files = os.listdir(path)
        score = False
        if os.path.lexists(path+os.sep+'fireball.in') and os.path.exists(path+os.sep+'answer.bas'):
            for ioutput in self.output_names:
                if os.path.isfile(path + os.sep + ioutput):
                    rf = open(path + os.sep + ioutput)
                    if "Welcome to FIREBALL" in rf.readline():
                        score = True
                    rf.close()

        if score:
            print(path)
            ret.append(path)

        for ifile in files:
            if os.path.isdir(path + os.sep + ifile):
                self.process_directory(path + os.sep + ifile, fireball_dirs)

    def save_list_candidates(self, filename):
        """
        Scan all databases looking for candidates for evaluation

        :return: A list of pairs, each pair contains the name of the database
                    and candidate identifier.
        """

        ret = []
        if not os.path.isfile(filename):
            self.process_directory(self.source_dir, ret)
            wf = open('fireball_directories.txt', 'w')
            for i in ret:
                wf.write("%s\n" % i)
            wf.close()
        else:
            rf = open(filename)
            ret = [x[:-1] for x in rf.readlines()]
            rf.close()

        return ret

    def is_evaluated(self, pcdb, path):
        if pcdb.entries.find({'properties.path': path}).count() > 0:
            return True
        else:
            return False

    def run(self):
        """
        Process the list of directories from 'source_dir' or 'source_file'
        Using 'nconcurrent' python processes.

        :return:
        """

        procs = []
        ids_running = []

        # Creating a list to store the 'nconcurrent' running jobs
        for i in range(self.nconcurrent):
            procs.append(None)
            ids_running.append(None)

        to_evaluate = []
        if self.source_file is None:
            self.process_directory(self.source_dir, to_evaluate)
        else:
            rf = open(self.source_file)
            to_evaluate = [x.strip() for x in rf.readlines()]

        # Main loop looking permanently for candidates for evaluation
        while True:

            index = 0
            currently_evaluating = 0
            for j in range(self.nconcurrent):
                if procs[j] is not None and procs[j].is_alive():
                    currently_evaluating += 1
            print('Candidates to evaluate: %d  Candidates in evaluation: %d' % (len(to_evaluate), currently_evaluating))

            while index < len(to_evaluate):

                pcdb = self.database

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
                    print('[Slot: %d] Evaluating: %s' % (slot, to_evaluate[index]))
                    procs[slot] = Process(target=self.worker, args=(to_evaluate[index],))
                    procs[slot].start()
                    time.sleep(0.5)
                else:
                    print('[Slot: %d] Evaluated: %s' % (slot, to_evaluate[index]))
                    pass

                index += 1
            time.sleep(self.sleeping_time)
