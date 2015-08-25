__author__ = 'Guillermo Avendano-Franco'

import os
import math
import socket
import time
import numpy as np
from multiprocessing import Process
from pychemia import pcm_log
from pychemia.analysis import StructureAnalysis
from pychemia.code.dftb.task import Relaxation
from pychemia.serializer import generic_serializer
from pychemia.db import get_database
from pychemia.utils.periodic import atomic_number
from pychemia.symm import StructureSymmetry


class EvaluatorDaemon:
    def __init__(self, database_settings, basedir, target_forces, nparal, relaxator_params,
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

    def run(self):

        def worker(db_settings, entry_id, workdir, target_forces, relaxator_params):

            pcdb = get_database(db_settings)

            pcm_log.info('[%s]: Starting relaxation. Target forces: %7.3e' % (str(entry_id), target_forces))

            if pcdb.is_locked(entry_id):
                return
            else:
                pcdb.lock(entry_id)
            structure = pcdb.get_structure(entry_id)
            structure = structure.scale()

            relaxer = Relaxation(structure, workdir, relaxator_params=relaxator_params,
                                        target_forces=target_forces, waiting=True)
            relaxer.run()
            pcm_log.info('[%s]: Finished relaxation. Target forces: %7.3e' % (str(entry_id), target_forces))

            filename = workdir + os.sep + 'detailed.out'
            if os.path.isfile(filename):
                forces, stress, total_energy = relaxer.get_forces_stress_energy()
                if forces is None:
                    pcm_log.error('No forces found on %s' % filename)
                if stress is None:
                    pcm_log.error('No stress found on %s' % filename)
                if total_energy is None:
                    pcm_log.error('No total_energy found on %s' % filename)

                new_structure = relaxer.get_final_geometry()

                if forces is not None and stress is not None and total_energy is not None and new_structure is not None:
                    pcm_log.info('[%s]: Updating properties' % str(entry_id))
                    symmetry = StructureSymmetry(new_structure)
                    pcdb.update(entry_id, structure=new_structure)
                    te = total_energy
                    pcdb.entries.update({'_id': entry_id},
                                        {'$set': {'status.relaxation': 'succeed',
                                                  'status.target_forces': target_forces,
                                                  'properties.forces': generic_serializer(forces),
                                                  'properties.stress': generic_serializer(stress),
                                                  'properties.energy': te,
                                                  'properties.energy_pa': te/new_structure.natom,
                                                  'properties.energy_pf': te/new_structure.get_composition().gcd,
                                                  'properties.spacegroup': symmetry.number()}})

                    # Fingerprint
                    # Update the fingerprints only if the two structures are really different
                    diffnatom = structure.natom != new_structure.natom
                    diffcell = np.max(np.abs((structure.cell - new_structure.cell).flatten()))
                    diffreduced = np.max(np.abs((structure.reduced - new_structure.reduced).flatten()))
                    if diffnatom != 0 or diffcell > 1E-7 or diffreduced > 1E-7:
                        analysis = StructureAnalysis(new_structure, radius=50)
                        x, ys = analysis.fp_oganov(delta=0.01, sigma=0.01)
                        fingerprint = {'_id': entry_id}
                        for k in ys:
                            atomic_number1 = atomic_number(new_structure.species[k[0]])
                            atomic_number2 = atomic_number(new_structure.species[k[1]])
                            pair = '%06d' % min(atomic_number1 * 1000 + atomic_number2,
                                                atomic_number2 * 1000 + atomic_number1)
                            fingerprint[pair] = list(ys[k])

                        if pcdb.db.fingerprints.find_one({'_id': entry_id}) is None:
                            pcdb.db.fingerprints.insert(fingerprint)
                        else:
                            pcdb.db.fingerprints.update({'_id': entry_id}, fingerprint)
                    else:
                        pcm_log.debug('Original and new structures are very similar.')
                        pcm_log.debug('Max diff cell: %10.3e' %
                                  np.max(np.absolute((structure.cell - new_structure.cell).flatten())))
                        if structure.natom == new_structure.natom:
                            pcm_log.debug('Max diff reduced coordinates: %10.3e' %
                                      np.max(np.absolute((structure.reduced - new_structure.reduced).flatten())))

                else:
                    pcdb.entries.update({'_id': entry_id}, {'$set': {'status.relaxation': 'failed'}})
                    pcm_log.error('Bad data after relaxation. Tagging relaxation as failed')
            else:
                pcm_log.error('ERROR: File not found %s' % filename)
            pcm_log.info('[%s]: Unlocking the entry' % str(entry_id))
            pcdb.unlock(entry_id)

        procs = []
        ids_running = []
        for i in range(self.nparal):
            procs.append(None)
            ids_running.append(None)

        pcdb = get_database(self.database_settings)
        for entry in pcdb.entries.find({}):
            pcdb.unlock(entry['_id'], name=socket.gethostname())
        print 'Number of entries: ', pcdb.entries.count()

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
            print '[%s]-> To Relax: %d  Relaxing: %d' % (pcdb.name, len(to_relax), sum(relaxing))
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
                            step_target = max([10 ** math.floor(math.log10(current_status / 2.0)), self.target_forces])
                            pcm_log.debug('New target  forces-stress: %7.3e' % step_target)

                            procs[j] = Process(target=worker, args=(db_settings, entry_id, workdir, step_target,
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

    @staticmethod
    def get_current_status(entry):
        if entry is not None and 'properties' in entry and entry['properties'] is not None:
            if 'forces' in entry['properties'] and entry['properties']['forces'] is not None:
                forces = np.max(np.abs(np.array(entry['properties']['forces'], dtype=float).flatten()))
            else:
                pcm_log.debug('No forces')
                forces = 10
            if 'stress' in entry['properties'] and entry['properties']['stress'] is not None:
                stress = np.max(np.abs(np.array(entry['properties']['stress'], dtype=float).flatten()))
            else:
                pcm_log.debug('No stress')
                stress = 10
        else:
            pcm_log.debug('Bad entry')
            print entry
            forces = 10
            stress = 10
        return max([forces, stress])

    def is_evaluable(self, entry):
        if 'lock' in entry['status']:
            return False
        elif self.evaluate_all:
            return True
        elif self.is_evaluated(entry):
            return False
        elif 'relaxation' not in entry['status']:
            return True
        if 'relaxation' in entry['status']:
            if entry['status']['relaxation'] == 'failed':
                if self.evaluate_failed:
                    return True
                else:
                    return False
        if self.get_current_status(entry) > self.target_forces:
            return True
        else:
            return False
