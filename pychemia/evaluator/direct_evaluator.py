from __future__ import print_function
import os
import socket
import time
from multiprocessing import Process
import numpy as np
from pychemia import pcm_log, HAS_PYMONGO

if HAS_PYMONGO:
    from pychemia.db import get_database


class DirectEvaluator:
    def __init__(self, database_settings, basedir, target_forces, nparal, relaxator_params, worker,
                 evaluate_failed=False, evaluate_all=False, waiting=False):

        self.database_settings = database_settings
        self.basedir = basedir
        self.nparal = nparal
        self.target_forces = target_forces
        self.relaxator_params = relaxator_params
        self.sleeping_time = 30
        self.evaluate_failed = evaluate_failed
        self.evaluate_all = evaluate_all
        self.waiting = waiting
        self.worker = worker

    def run(self):

        procs = []
        ids_running = []
        for i in range(self.nparal):
            procs.append(None)
            ids_running.append(None)

        pcdb = get_database(self.database_settings)
        print('Database contains: %d entries' % pcdb.entries.count())

        for entry in pcdb.entries.find({}):
            pcdb.unlock(entry['_id'], name=socket.gethostname())
        print('Number of entries: ', pcdb.entries.count())

        while True:

            to_relax = []

            for entry in pcdb.entries.find({}):
                if self.is_evaluable(entry):
                    pcm_log.debug('Adding entry %s from db %s' % (str(entry['_id']), pcdb.name))
                    to_relax.append(entry['_id'])

            if len(to_relax) == 0:
                pcm_log.debug('No more entries to evaluate, sleeping %d seconds' % self.sleeping_time)
                time.sleep(self.sleeping_time)
            else:
                sortlist = np.array([pcdb.get_structure(i).natom for i in to_relax]).argsort()[::-1]
                to_relax = list(np.array(to_relax)[sortlist])
                pcm_log.debug('Number of entries to evaluate: %d' % len(to_relax))

            index = 0
            relaxing = np.array([True for j in range(self.nparal) if procs[j] is not None and procs[j].is_alive()])
            print('[%s]-> To Relax: %d  Relaxing: %d' % (pcdb.name, len(to_relax), sum(relaxing)))
            while index < len(to_relax):

                entry_id = to_relax[index]
                entry = pcdb.entries.find_one({'_id': entry_id})
                if not os.path.exists(self.basedir + os.sep + pcdb.name):
                    os.mkdir(self.basedir + os.sep + pcdb.name)
                for j in range(self.nparal):
                    if procs[j] is None or not procs[j].is_alive():
                        if ids_running[j] is not None:
                            pcm_log.debug('%s is not alive. Exit code: %2d. Locked: %s' %
                                          (str(ids_running[j]), procs[j].exitcode, str(pcdb.is_locked(ids_running[j]))))
                        if self.is_evaluable(pcdb.entries.find_one({'_id': entry_id})):
                            pcm_log.debug('Evaluable: %s. Relaxing entry %d of %d' %
                                          (str(entry_id), index, len(to_relax)))
                            ids_running[j] = entry_id
                            workdir = self.basedir + os.sep + pcdb.name + os.sep + str(entry_id)
                            if not os.path.exists(workdir):
                                os.mkdir(workdir)
                            db_settings = self.database_settings.copy()
                            pcm_log.debug('Launching for %s id: %s' % (pcdb.name, str(entry_id)))

                            # Relax lowering the target forces by one order of magnitude each time
                            current_status = self.get_current_status(entry)
                            pcm_log.debug('Current max forces-stress: %7.3e' % current_status)
                            step_target = max(0.1 * current_status, self.target_forces)
                            pcm_log.debug('New target  forces-stress: %7.3e' % step_target)

                            procs[j] = Process(target=self.worker, args=(db_settings, entry_id, workdir, step_target,
                                                                         self.relaxator_params))
                            procs[j].start()
                        else:
                            pcm_log.debug('Not evaluable: %s' % str(entry_id))

                        index += 1
                        break
                time.sleep(3)
            time.sleep(self.sleeping_time)

    def is_evaluated(self, entry):
        return self.get_current_status(entry) < self.target_forces

    def get_current_status(self, entry):
        if entry is not None and 'properties' in entry and entry['properties'] is not None:
            if 'forces' in entry['properties'] and entry['properties']['forces'] is not None:
                forces = np.array(entry['properties']['forces']).reshape((-1, 3))
                max_force = np.max(np.apply_along_axis(np.linalg.norm, 1, forces))
            else:
                pcm_log.debug('No forces')
                max_force = 1
            if 'stress' in entry['properties'] and entry['properties']['stress'] is not None:
                stress = np.array(entry['properties']['stress']).reshape((-1, 3))
                max_stress = np.max(np.abs(stress.flatten()))
            else:
                pcm_log.debug('No stress')
                max_stress = 1
        else:
            pcm_log.debug('Bad entry')
            print(entry)
            max_force = 1
            max_stress = 1
        if max(max_force, max_stress) > self.target_forces:
            print('Status for %s: forces: %9.2E stress: %9.2E' % (entry['_id'], max_force, max_stress))
        return max(max_force, max_stress)

    def is_evaluable(self, entry):
        if 'lock' in entry['status']:
            return False
        elif self.evaluate_all:
            return True
        elif self.get_current_status(entry) > self.target_forces:
            return True
        else:
            return False
