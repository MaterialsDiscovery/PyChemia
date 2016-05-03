#!/usr/bin/env python

import argparse
import logging
import os

import numpy as np

from pychemia import pcm_log
from pychemia.analysis import StructureAnalysis
from pychemia.code.dftb import read_detailed_out
from pychemia.code.dftb.task import Relaxation
from pychemia.db import get_database
from pychemia.evaluator import DirectEvaluator
from pychemia.serializer import generic_serializer
from pychemia.utils.periodic import atomic_number


def worker(db_settings, entry_id, workdir, target_forces, relaxator_params):
    pcdb = get_database(db_settings)
    pcm_log.info('[%s]: Starting relaxation. Target forces: %7.3e' % (str(entry_id), target_forces))

    if pcdb.is_locked(entry_id):
        return
    else:
        pcdb.lock(entry_id)
    structure = pcdb.get_structure(entry_id)
    structure = structure.scale()
    if 'forced' in relaxator_params:
        forced = relaxator_params['forced']
    else:
        forced = True

    relaxer = Relaxation(structure, relaxator_params=relaxator_params, workdir=workdir,
                         target_forces=target_forces, waiting=True, forced=forced)
    relaxer.run()
    pcm_log.info('[%s]: Finished relaxation. Target forces: %7.3e' % (str(entry_id), target_forces))

    filename = workdir + os.sep + 'detailed.out'
    if os.path.isfile(filename):
        detailed = read_detailed_out(filename)
        if 'max_force' in detailed:
            max_force = detailed['max_force']
        else:
            max_force = target_forces

        if 'max_deriv' in detailed:
            max_deriv = detailed['max_deriv']
        else:
            max_deriv = target_forces

        forces, stress, total_energy = relaxer.get_forces_stress_energy()

        if forces is not None:
            magnitude_forces = np.apply_along_axis(np.linalg.norm, 1, forces)
            print('Forces: Max: %9.3e Avg: %9.3e' % (np.max(magnitude_forces), np.average(magnitude_forces)))
            print('Stress: ', np.max(np.abs(stress.flatten())))

        if forces is None:
            pcm_log.error('No forces found on %s' % filename)
        if stress is None:
            pcm_log.error('No stress found on %s' % filename)
        if total_energy is None:
            pcm_log.error('No total_energy found on %s' % filename)

        new_structure = relaxer.get_final_geometry()

        if forces is not None and stress is not None and total_energy is not None and new_structure is not None:
            pcm_log.info('[%s]: Updating properties' % str(entry_id))
            pcdb.update(entry_id, structure=new_structure)
            te = total_energy
            pcdb.entries.update({'_id': entry_id},
                                {'$set': {'status.relaxation': 'succeed',
                                          'status.target_forces': target_forces,
                                          'properties.forces': generic_serializer(forces),
                                          'properties.stress': generic_serializer(stress),
                                          'properties.energy': te,
                                          'properties.max_force': max_force,
                                          'properties.max_deriv': max_deriv,
                                          'properties.energy_pa': te / new_structure.natom,
                                          'properties.energy_pf': te / new_structure.get_composition().gcd}})

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
                pcm_log.debug('Max diff cell: %10.3e' % np.max(np.absolute((structure.cell -
                                                                            new_structure.cell).flatten())))
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


if __name__ == '__main__':

    logging.basicConfig(level=logging.DEBUG)
    logger = logging.getLogger('pychemia')
    logger.addHandler(logging.NullHandler())
    logger.setLevel(logging.DEBUG)

    description = """Launch DFTB+ for non-evaluated entries in a PyChemia Database"""

    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-t', '--host',
                        default='localhost', metavar='server', type=str,
                        help='Hostname or address (default: localhost)')
    parser.add_argument('-o', '--port',
                        default=27017, metavar='port', type=int,
                        help='MongoDB port (default: 27017)')
    parser.add_argument('-u', '--user',
                        default=None, metavar='username', type=str,
                        help='Username (default: None)')
    parser.add_argument('-p', '--passwd',
                        default=None, metavar='password', type=str,
                        help='Password (default: None)')
    parser.add_argument('-d', '--dbname',
                        default=None, metavar='dbname', type=str,
                        help='PyChemia Database name (default: None)')
    parser.add_argument('-l', '--slater_path',
                        default=None, metavar='path', type=str,
                        help='Slater Path (default: None)')
    parser.add_argument('-f', '--target_forces',
                        default=1E-3, metavar='x', type=float,
                        help='Target Forces (default: 1E-3)')
    parser.add_argument('-n', '--nparal',
                        default=1, metavar='N', type=int,
                        help='Number of parallel processes (default: 1)')
    parser.add_argument('-r', '--replicaset',
                        default=None, metavar='name', type=str,
                        help='ReplicaSet  (default: None)')
    parser.add_argument('-w', '--workdir',
                        default='.', metavar='path', type=str,
                        help='Working Directory  (default: None)')
    parser.add_argument('--evaluate_all', action='store_true',
                        help='Evaluate All  (default: No)')
    parser.add_argument('--waiting', action='store_true',
                        help='Waiting  (default: No)')
    parser.add_argument('--ssl', action='store_true',
                        help='Use SSL to connect to MongoDB  (default: No)')

    args = parser.parse_args()
    if args.dbname is None:
        parser.print_help()
        exit(1)

    print(args)
    db_settings = {'name': args.dbname, 'host': args.host, 'port': args.port, 'ssl': args.ssl,
                   'replicaset': args.replicaset}
    if args.user is not None:
        if args.passwd is None:
            raise ValueError('Password is mandatory if user is entered')
        db_settings['user'] = args.user
        db_settings['passwd'] = args.passwd
    relaxator_params = {'slater_path': args.slater_path}
    print('pyChemia Evaluator using DFTB+')
    print('dbname    : %s' % args.dbname)
    print('host      : %s' % args.host)
    print('port      : %d' % args.port)
    print('user      : %s' % args.user)
    print('replicaset: %s' % args.replicaset)
    print('workdir   : %s' % args.workdir)
    print('nparal    : %d' % args.nparal)
    print('slater-path   : %s' % str(args.slater_path))
    print('target-forces : %.2E' % args.target_forces)
    print('evaluate_all  : %s' % str(args.evaluate_all))
    print('waiting       : %s' % str(args.waiting))
    print('ssl           : %s' % str(args.ssl))

    print(db_settings)
    evaluator = DirectEvaluator(db_settings, args.workdir,
                                target_forces=args.target_forces,
                                nparal=args.nparal,
                                relaxator_params=relaxator_params,
                                worker=worker,
                                evaluate_all=args.evaluate_all,
                                waiting=args.waiting)
    evaluator.run()
