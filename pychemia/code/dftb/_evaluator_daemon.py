__author__ = 'Guillermo Avendano-Franco'

import os
import math
import socket
import time
import numpy as np
from multiprocessing import Process
from pychemia import log
from pychemia.analysis import StructureAnalysis
from pychemia.code.dftb import Relaxator
from pychemia.serializer import generic_serializer
from pychemia.db import get_database
from pychemia.utils.periodic import atomic_number
from pychemia.symm import StructureSymmetry


class EvaluatorDaemon():
    def __init__(self, database_settings, basedir, target_forces, nparal, slater_path,
                 evaluate_failed=False, evaluate_all=False):

        self.database_settings = database_settings
        self.basedir = basedir
        self.nparal = nparal
        self.target_forces = target_forces
        if isinstance(slater_path, basestring):
            self.slater_path = [slater_path]
        for i in self.slater_path:
            assert os.path.isdir(i)
        self.sleeping_time = 30
        self.evaluate_failed = evaluate_failed
        self.evaluate_all = evaluate_all

    def run(self):

        def worker(db_settings, entry_id, workdir, target_forces, slater_path):

            db = get_database(db_settings)

            log.info('Starting relaxation for %s with target forces of %7.3e' % (str(entry_id), target_forces))

            if db.is_locked(entry_id):
                return
            else:
                db.lock(entry_id)
            structure = db.get_structure(entry_id)
            structure_dict, properties, status = db.get_dicts(entry_id)

            if properties is None:
                properties = {}
            if status is None:
                status = {}

            relaxer = Relaxator(workdir, structure, slater_path, target_forces, waiting=True)
            relaxer.run()
            log.info('Finished relaxation for %s with target forces of %7.3e' % (str(entry_id), target_forces))

            filename = workdir + os.sep + 'detailed.out'
            if os.path.isfile(filename):
                forces, stress, total_energy = relaxer.get_forces_stress_energy()
                if forces is None:
                    log.error('No forces found on %s' % filename)
                if stress is None:
                    log.error('No stress found on %s' % filename)
                if total_energy is None:
                    log.error('No total_energy found on %s' % filename)

                new_structure = relaxer.get_final_geometry()

                if forces is not None and stress is not None and total_energy is not None and new_structure is not None:
                    log.info('Updating the database with new properties: %s' % str(entry_id))
                    status['relaxation'] = 'succeed'
                    status['target_forces'] = target_forces
                    properties['forces'] = generic_serializer(forces)
                    properties['stress'] = generic_serializer(stress),
                    properties['energy'] = total_energy
                    properties['energy_pa'] = total_energy / new_structure.natom
                    properties['energy_pf'] = total_energy / new_structure.get_composition().gcd
                    symmetry = StructureSymmetry(new_structure)
                    properties['spacegroup'] = symmetry.number()
                    db.update(entry_id, structure=new_structure, properties=properties, status=status)

                    # Fingerprint
                    # Update the fingerprints only if the two structures are really different
                    if structure.natom != new_structure.natom or \
                                    np.max(np.absolute((structure.cell - new_structure.cell).flatten())) > 1E-7 or \
                                    np.max(np.absolute((structure.reduced - new_structure.reduced).flatten())) > 1E-7:

                        analysis = StructureAnalysis(new_structure, radius=50)
                        x, ys = analysis.fp_oganov(delta=0.01, sigma=0.01)
                        fingerprint = {'_id': entry_id}
                        for k in ys:
                            atomic_number1 = atomic_number(new_structure.species[k[0]])
                            atomic_number2 = atomic_number(new_structure.species[k[1]])
                            pair = '%06d' % min(atomic_number1 * 1000 + atomic_number2,
                                                atomic_number2 * 1000 + atomic_number1)
                            fingerprint[pair] = list(ys[k])

                        if db.db.fingerprints.find_one({'_id': entry_id}) is None:
                            db.db.fingerprints.insert(fingerprint)
                        else:
                            db.db.fingerprints.update({'_id': entry_id}, fingerprint)
                    else:
                        log.debug('Original and new structures are very similar.')
                        log.debug('Max diff cell: %10.3e' %
                                  np.max(np.absolute((structure.cell - new_structure.cell).flatten())))
                        if structure.natom == new_structure.natom:
                            log.debug('Max diff reduced coordinates: %10.3e' %
                                      np.max(np.absolute((structure.reduced - new_structure.reduced).flatten())))

                else:
                    status['relaxation'] = 'failed'
                    db.update(entry_id, status=status)
                    log.error('Bad data after relaxation. Tagging relaxation as failed')
            else:
                log.error('ERROR: File not found %s' % filename)
            log.info('Unlocking the entry: %s' % str(entry_id))
            db.unlock(entry_id)

        procs = []
        ids_running = []
        for i in range(self.nparal):
            procs.append(None)
            ids_running.append(None)

        db = get_database(self.database_settings)
        for entry in db.entries.find({}):
            db.unlock(entry['_id'], name=socket.gethostname())

        while True:

            to_relax = []

            for entry in db.entries.find({}):
                if self.is_evaluable(entry):
                    log.debug('Adding entry %s from db %s' % (str(entry['_id']), db.name))
                    to_relax.append(entry['_id'])

            if len(to_relax) == 0:
                log.debug('No more entries to evaluate, sleeping %d seconds' % self.sleeping_time)
                time.sleep(self.sleeping_time)
            else:
                sortlist = np.array([db.get_structure(i).natom for i in to_relax]).argsort()[::-1]
                to_relax = list(np.array(to_relax)[sortlist])
                log.debug('Number of entries to evaluate: %d' % len(to_relax))

            index = 0
            while index < len(to_relax):

                entry_id = to_relax[index]
                entry = db.entries.find_one({'_id': entry_id})
                if not os.path.exists(self.basedir + os.sep + db.name):
                    os.mkdir(self.basedir + os.sep + db.name)
                for j in range(self.nparal):
                    if procs[j] is None or not procs[j].is_alive():
                        if ids_running[j] is not None:
                            log.debug('%s is not alive. Exit code: %2d. Locked: %s' %
                                      (str(ids_running[j]), procs[j].exitcode, str(db.is_locked(ids_running[j]))))
                        if self.is_evaluable(db.entries.find_one({'_id': entry_id})):
                            log.debug('Evaluable: %s' % str(entry_id))
                            ids_running[j] = entry_id
                            workdir = self.basedir + os.sep + db.name + os.sep + str(entry_id)
                            if not os.path.exists(workdir):
                                os.mkdir(workdir)
                            db_settings = self.database_settings.copy()
                            log.debug('Launching for %s id: %s' % (db.name, str(entry_id)))

                            # Relax lowering the target forces by one order of magnitude each time
                            current_status = self.get_current_status(entry)
                            log.debug('Current max forces-stress: %7.3e' % current_status)
                            step_target = max([10 ** math.floor(math.log10(current_status / 2.0)), self.target_forces])
                            log.debug('New target  forces-stress: %7.3e' % step_target)

                            procs[j] = Process(target=worker, args=(db_settings, entry_id, workdir, step_target,
                                                                    self.slater_path))
                            procs[j].start()
                        else:
                            log.debug('Not evaluable: %s' % str(entry_id))

                        index += 1
                        break
                time.sleep(1)
            time.sleep(self.sleeping_time)

    def is_evaluated(self, entry):
        return self.get_current_status(entry) < self.target_forces

    def get_current_status(self, entry):
        if entry is not None and 'properties' in entry and entry['properties'] is not None:
            if 'forces' in entry['properties'] and entry['properties']['forces'] is not None:
                forces = np.max(np.abs(np.array(entry['properties']['forces'], dtype=float).flatten()))
            else:
                log.debug('No forces')
                forces = 10
            if 'stress' in entry['properties'] and entry['properties']['stress'] is not None:
                stress = np.max(np.abs(np.array(entry['properties']['stress'], dtype=float).flatten()))
            else:
                log.debug('No stress')
                stress = 10
        else:
            log.debug('Bad entry')
            print entry
            forces = 10
            stress = 10
        return max([forces, stress])

    def is_evaluable(self, entry):
        if 'status' in entry and 'lock' in entry['status']:
            return False
        elif self.evaluate_all:
            return True
        elif self.is_evaluated(entry):
            return False
        elif 'status' not in entry or entry['status'] is None or 'relaxation' not in entry['status']:
            return True
        if 'status' in entry and entry['status'] is not None and 'relaxation' in entry['status']:
            if entry['status']['relaxation'] == 'failed' and self.evaluate_failed:
                return True
        if self.get_current_status(entry) > self.target_forces:
            return True
        else:
            return False
