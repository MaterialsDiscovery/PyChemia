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
                 evaluate_failed=False, evaluate_all=False, waiting=False, pressure=0.0):

        self.database_settings = database_settings
        self.basedir = basedir
        self.nparal = nparal
        self.target_forces = target_forces
        self.relaxator_params = relaxator_params
        self.sleeping_time = 120
        self.evaluate_failed = evaluate_failed
        self.evaluate_all = evaluate_all
        self.waiting = waiting
        self.worker = worker
        # Pressure in kB
        self.pressure = pressure

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

    def get_current_status(self, entry, verbose=False):

        max_force = 1
        max_diag_stress = 1
        max_nondiag_stress = 1
        
        if entry is not None and 'properties' in entry and entry['properties'] is not None:
            if 'forces' in entry['properties'] and entry['properties']['forces'] is not None:
                forces = np.array(entry['properties']['forces']).reshape((-1, 3))
                max_force = np.max(np.apply_along_axis(np.linalg.norm, 1, forces))
            if 'stress' in entry['properties'] and entry['properties']['stress'] is not None:
                stress = np.array(entry['properties']['stress']).reshape((-1, 3))
                assert(stress.shape==(2,3))
                max_diag_stress=np.abs(np.max(np.abs(stress[0])) - (self.pressure/1602.1766208))
                max_nondiag_stress=np.max(np.abs(stress[1]))
        else:
            pcm_log.debug('Bad entry')
            print(entry)
        if max(max_force, max_diag_stress, max_nondiag_stress) > self.target_forces and verbose:
            if (max_force*max_diag_stress*max_nondiag_stress) == 1.0:
                print('No forces/stress information for entry: %s' % entry['_id'])
            else:
                print('Convergence status for entry: %s' % entry['_id'])
                print('Max Interatomic Force: %9.2E [eV/Ang]' % max_force)
                print('Max Diag stress      : %9.2E with pressure: %9.2E [eV/Ang^3]' % (max_diag_stress, self.pressure/1602.1766208 ))
                print('Max NonDiag stress   : %9.2E [eV/Ang^3]' % max_nondiag_stress)
        return max(max_force, max_diag_stress, max_nondiag_stress)

    def is_evaluable(self, entry):
        if 'lock' in entry['status']:
            print('Entry is not evaluable because %s is locked' % entry['_id'])
            return False
        elif self.evaluate_all:
            print('Entry is evaluable because evaluate all')
            return True
        elif self.get_current_status(entry, verbose=True) > self.target_forces:
            print('Entry is evaluable because forces not converged %f' % self.get_current_status(entry))
            return True
        else:
            return False
