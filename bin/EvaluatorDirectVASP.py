#!/usr/bin/env python

import argparse
import logging
import os
import re
import shutil
import numpy as np
import subprocess
from pychemia import pcm_log, Structure
from pychemia.analysis import StructureAnalysis
from pychemia.code.vasp import write_poscar, read_poscar
from pychemia.code.vasp.task import IonRelaxation
from pychemia.db import get_database
from pychemia.evaluator import DirectEvaluator
from pychemia.utils.serializer import generic_serializer
from pychemia.utils.periodic import atomic_number


def worker_maise(db_settings, entry_id, workdir, target_forces, relaxator_params):

    max_ncalls = 6
    pcdb = get_database(db_settings)
    pcm_log.info('[%s]: Starting relaxation. Target forces: %7.3e' % (str(entry_id), target_forces))

    if pcdb.is_locked(entry_id):
        return
    else:
        pcdb.lock(entry_id)
    structure = pcdb.get_structure(entry_id)
    status = pcdb.get_dicts(entry_id)[2]

    if 'ncalls' in status:
        ncalls = status['ncalls'] + 1 
    else:
        ncalls = 1 

    # print('Current directory: '+os.getcwd() )
    # print('Working directory: '+workdir)
    write_poscar(structure, workdir+os.sep+'POSCAR')
    if not os.path.exists(workdir+os.sep+'setup'):
        shutil.copy2('setup', workdir)
    if not os.path.exists(workdir+os.sep+'INI'):
        os.symlink(os.getcwd()+os.sep+'INI', workdir+os.sep+'INI')
    if not os.path.exists(workdir+os.sep+'maise'):
        os.symlink(os.getcwd()+os.sep+'maise', workdir+os.sep+'maise')
    cwd = os.getcwd()
    os.chdir(workdir)
    wf = open('maise.stdout', 'w')
    subprocess.call(['./maise'], stdout=wf)
    wf.close()
    if os.path.isfile('OSZICAR'):
        energies = np.loadtxt('OSZICAR')
    else:
        energies = None
    if os.path.isfile('OUTCAR'):
        rf = open('OUTCAR', 'r')
        data = rf.read()
        state_dist = re.findall(r'Total CPU', data)                                       # WIH
        
        pos_forces = re.findall(r'TOTAL-FORCE \(eV/Angst\)\s*-*\s*([-.\d\s]+)\s+-{2}', data)
        pos_forces = np.array([x.split() for x in pos_forces], dtype=float)

        if len(pos_forces) > 0 and len(pos_forces[-1]) % 7 == 0:
            pos_forces.shape = (len(pos_forces), -1, 7)
            forces = pos_forces[:, :, 3:6]
            positions = pos_forces[:, :, :3]
        else:
            print('Forces and Positions could not be parsed : ', pos_forces.shape)
            print('pos_forces =\n%s ' % pos_forces)
            
        str_stress = re.findall('Total([\.\d\s-]*)in', data)
        if len(str_stress) == 2:
            stress = np.array([[float(y) for y in x.split()] for x in str_stress])
        else:
            stress = None
        str_stress = re.findall('in kB([\.\d\s-]*)energy', data)
        if len(str_stress) == 2:
            stress_kb = np.array([[float(y) for y in x.split()] for x in str_stress])
        else:
            stress_kb = None
    else:
        forces = None
        stress = None
        stress_kb = None
    state_dist = 'F'                                                                      # WIH

    new_structure = read_poscar('CONTCAR')
    if len(state_dist) == 0:                                                              # WIH
        print('WARNING: MAISE found distances too short in structure:', entry_id)           # WIH
        new_structure = Structure.random_cell(structure.composition)                          # WIH
    #if np.min(new_structure.distance_matrix()+np.eye(new_structure.natom))<0.23:       # WIH
    #    print('WARNING: Structure collapse 2 atoms, creating a new random structure')  # WIH
    #    new_structure=Structure.random_cell(new_structure.composition)                 # WIH
    if ncalls > max_ncalls:
        print('WARNING: Too many calls to MAISE and no relaxation succeeded, replacing structure: ', entry_id)    # WIH
        new_structure = Structure.random_cell(structure.composition)                      #WIH
        pcdb.entries.update({'_id': entry_id}, {'$set': {'status.ncalls': 0}})
    else:
        pcdb.entries.update({'_id': entry_id}, {'$set': {'status.ncalls': ncalls}})
    pcdb.update(entry_id, structure=new_structure)

    if energies is not None and forces is not None and stress is not None:

        te = energies[1]
        pcdb.entries.update({'_id': entry_id},
                            {'$set': {'status.relaxation': 'succeed',
                                      'status.target_forces': target_forces,
                                      'properties.initial_forces': generic_serializer(forces[0]),
                                      'properties.initial_stress': generic_serializer(stress[0]),
                                      'properties.initial_stress_kB': generic_serializer(stress_kb[0]),
                                      'properties.forces': generic_serializer(forces[1]),
                                      'properties.stress': generic_serializer(stress[1]),
                                      'properties.stress_kB': generic_serializer(stress_kb[1]),
                                      'properties.energy': te,
                                      'properties.energy_pa': te / new_structure.natom,
                                      'properties.energy_pf': te / new_structure.get_composition().gcd}})

    for ifile in ['POSCAR', 'CONTCAR', 'setup', 'OUTCAR', 'maise.stdout', 'list.dat']:
        if not os.path.exists(ifile):
            wf = open(ifile, 'w')
            wf.write('')
            wf.close()
        n = 1
        while True:
            if os.path.exists(ifile + ('_%03d' % n)):
                n += 1
            else:
                break
        os.rename(ifile, ifile+('_%03d' % n))

    pcm_log.info('[%s]: Unlocking the entry' % str(entry_id))
    pcdb.unlock(entry_id)


def worker_vasp(db_settings, entry_id, workdir, target_forces, relaxator_params):
    pcdb = get_database(db_settings)
    pcm_log.info('[%s]: Starting relaxation. Target forces: %7.3e' % (str(entry_id), target_forces))

    if pcdb.is_locked(entry_id):
        return
    else:
        pcdb.lock(entry_id)
    structure = pcdb.get_structure(entry_id)
    structure = structure.scale()
    print('relaxator_params', relaxator_params)
    relaxer = IonRelaxation(structure, workdir=workdir, target_forces=target_forces, waiting=False,
                            binary=relaxator_params['binary'], encut=1.3, kp_grid=None, kp_density=1E4,
                            relax_cell=True, max_calls=10)
    print('relaxing on:', relaxer.workdir)
    relaxer.run(relaxator_params['nmpiparal'])
    pcm_log.info('[%s]: Finished relaxation. Target forces: %7.3e' % (str(entry_id), target_forces))

    filename = workdir + os.sep + 'OUTCAR'
    if os.path.isfile(filename):

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

    description = """Launch VASP for non-evaluated entries in a PyChemia Database"""

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
    parser.add_argument('-b', '--binary',
                        default=None, metavar='path', type=str,
                        help='VASP binary (default: None)')
    parser.add_argument('-f', '--target_forces',
                        default=1E-3, metavar='x', type=float,
                        help='Target Forces (default: 1E-3)')
    parser.add_argument('--pressure_kB',
                        default=0.0, metavar='x', type=float,
                        help='Imposed Pressure (default: 0.0 kB)')
    parser.add_argument('-n', '--nparal',
                        default=1, metavar='N', type=int,
                        help='Number of parallel processes (default: 1)')
    parser.add_argument('-m', '--nmpiparal',
                        default=1, metavar='N', type=int,
                        help='Number of MPI parallel processes (default: 1)')
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
    if args.binary is None:
        args.binary = 'vasp'
    if args.nmpiparal is None:
        args.nmpiparal = 1
    relaxator_params = {'binary': args.binary, 'nmpiparal': args.nmpiparal}
    print('pyChemia Evaluator using VASP')
    print('dbname    : %s' % args.dbname)
    print('host      : %s' % args.host)
    print('port      : %d' % args.port)
    print('user      : %s' % args.user)
    print('replicaset: %s' % args.replicaset)
    print('workdir   : %s' % args.workdir)
    print('nparal    : %d' % args.nparal)
    print('nmpiparal : %d' % args.nmpiparal)
    print('binary    : %s' % str(args.binary))
    print('target-forces : %.2E' % args.target_forces)
    print('pressure_kB   : %.2E' % args.pressure_kB)
    print('evaluate_all  : %s' % str(args.evaluate_all))
    print('waiting       : %s' % str(args.waiting))
    print('ssl           : %s' % str(args.ssl))

    print(db_settings)

    if args.binary[:4].lower() == 'vasp':

        evaluator = DirectEvaluator(db_settings, args.workdir,
                                    target_forces=args.target_forces,
                                    nparal=args.nparal,
                                    relaxator_params=relaxator_params,
                                    worker=worker_vasp,
                                    evaluate_all=args.evaluate_all,
                                    waiting=args.waiting,
                                    pressure=args.pressure_kB)
    elif args.binary[:4].lower() == 'mais':
        evaluator = DirectEvaluator(db_settings, args.workdir,
                                    target_forces=args.target_forces,
                                    nparal=args.nparal,
                                    relaxator_params=relaxator_params,
                                    worker=worker_maise,
                                    evaluate_all=args.evaluate_all,
                                    waiting=args.waiting,
                                    pressure=args.pressure_kB)
    else:
        print('I could not recognize the binary you pretent to use')
        print('Use a binary such as vasp* or maise')
        exit(1)

    evaluator.run()
