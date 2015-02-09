__author__ = 'Guillermo Avendano-Franco'

import os
from pychemia import Structure,log
from pychemia.code.dftb import Relaxator, read_geometry_gen, read_detailed_out
from pychemia.serializer import generic_serializer
from pychemia.db import PyChemiaDB
from multiprocessing import Process
import time
import numpy as np

class EvaluatorDaemon():

    def __init__(self, databases, basedir, target_forces, nparal, slater_path):

        self.databases = databases
        self.basedir = basedir
        self.nparal = nparal
        self.target_forces = target_forces
        self.slater_path = slater_path

    def run(self):

        def worker(dbname, entry_id, workdir, target_forces, slater_path):

            db = PyChemiaDB(dbname)
            entry = db.entries.find_one({'_id': entry_id})
            structure = Structure.from_dict(entry['structure'])
            relaxer=Relaxator(workdir, structure, slater_path, target_forces)
            relaxer.run()

            filename = workdir+os.sep+'detailed.out'
            if os.path.isfile(filename):
                forces, stress, total_energy = read_detailed_out(filename=filename)
                geometry = workdir+os.sep+'geo_end.gen'
                if forces is not None and stress is not None and total_energy is not None and os.path.isfile(geometry):
                    new_structure = read_geometry_gen(geometry)
                    properties = {'forces': generic_serializer(forces), 'stress': generic_serializer(stress), 'energy': total_energy}
                    entry['structure'] = new_structure.to_dict()
                    entry['properties'] = properties
                    db.update(entry_id, entry)

        procs = []
        for i in range(self.nparal):
            procs.append(None)

        while True:

            to_relax = []

            for idb in self.databases:

                db = PyChemiaDB(idb)
                for entry in db.entries.find({}):
                    if not self.is_evaluated(entry):
                        entry_id = entry['_id']
                        log.debug('Adding entry %s from db %s' % (str(entry_id), idb))
                        to_relax.append((idb, entry_id))

            if len(to_relax) == 0:
                log.debug('No more entries to evaluate, exiting')
                break

            index = 0
            while index < len(to_relax):

                idb = to_relax[index][0]
                entry_id = to_relax[index][1]
                if not os.path.exists(self.basedir + os.sep + idb):
                    os.mkdir(self.basedir + os.sep + idb)
                for j in range(self.nparal):
                    if procs[j] is None or not procs[j].is_alive():
                        log.debug('Relaxing for %s id: %s' % (idb, str(entry_id)))
                        workdir = self.basedir + os.sep + idb + os.sep + str(entry_id)
                        if not os.path.exists(workdir):
                            os.mkdir(workdir)
                        procs[j] = Process(target=worker, args=(idb, entry_id, workdir, self.target_forces, self.slater_path))
                        procs[j].start()
                        index += 1
                        break
                time.sleep(2)

    def is_evaluated(self, entry):
        if entry is not None and entry['properties'] is not None:
            if 'energy' not in entry['properties']:
                return False
            if 'stress' not in entry['properties']:
                return False
            if 'forces' not in entry['properties']:
                return False
            if 'forces' in entry['properties'] and np.max(np.abs(np.array(entry['properties']['forces'], dtype=float).flatten())) > self.target_forces:
                return False
            if 'stress' in entry['properties'] and np.max(np.abs(np.array(entry['properties']['stress'], dtype=float).flatten())) > self.target_forces:
                return False
            return True
        else:
            return False
