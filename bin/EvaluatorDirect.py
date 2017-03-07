#!/usr/bin/env python

import os
import re
import shutil
import argparse
import logging
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
from pychemia.code.vasp.task import IonRelaxation2
from pychemia.code.dftb.task import Relaxation
from pychemia.code.dftb import read_detailed_out


def worker_maise(db_settings, entry_id, workdir, relaxator_params):
    """
    Relax and return evaluate the energy of the structure stored with identifier 'entry_id'
     using the MAISE code

    :param db_settings: (dict) Dictionary of DB parameters needed to create a PyChemiaDB object
    :param entry_id: MongoDB identifier of one entry of the database created from db_settings
    :param workdir: (str) Working directory where input and output from MAISE is written
    :param relaxator_params: (dict) Arguments needed to control the relaxation using MAISE
                                Arguments are store as keys and they include:
                                'target_forces' : Used to defined the tolerance to consider one candidate as relaxed.
                                'source_dir': Directory with executable maise and directory INI
    :return:
    """
    max_ncalls = 6
    pcdb = get_database(db_settings)
    target_forces = relaxator_params['target_forces']
    source_dir = relaxator_params['source_dir']

    pcm_log.info('[%s]: Starting relaxation. Target forces: %7.3e' % (str(entry_id), target_forces))

    if pcdb.is_locked(entry_id):
        return
    else:
        pcdb.lock(entry_id)
    structure = pcdb.get_structure(entry_id)
    status = pcdb.get_dicts(entry_id)[2]

    if 'ncalls' in status and status['ncalls'] > 0:
        ncalls = status['ncalls'] + 1
        print('ncalls = ', status['ncalls'])
    else:
        ncalls = 1
    print('Verifing initial structure...')
    while np.min(structure.distance_matrix()+(np.eye(structure.natom)*5)) < 1.9:
        print('ERROR: Bad initial guess, two atoms are to close. Creating new random structure for id: %s' %
              str(entry_id))
        write_poscar(structure, workdir + os.sep + 'Fail_initial_POSCAR')  # WIH
        structure = Structure.random_cell(structure.composition)

    write_poscar(structure, workdir + os.sep + 'POSCAR')
    if not os.path.exists(workdir + os.sep + 'setup') and ncalls == 1:     # WIH
        print('First run.')  # WIH
        #   print('Verifying that everything runs smoothly') # WIH
        print(workdir + os.sep + 'setup')
        shutil.copy2(source_dir + os.sep + 'setup_1', workdir + os.sep + 'setup')   # WIH
    elif ncalls > 1:  # WIH
        shutil.copy2(source_dir + os.sep + 'setup_2', workdir + os.sep + 'setup')   # WIH
    if not os.path.exists(workdir + os.sep + 'INI'):
        os.symlink(source_dir + os.sep + 'INI', workdir + os.sep + 'INI')
    if not os.path.exists(workdir + os.sep + 'maise'):
        os.symlink(source_dir + os.sep + 'maise', workdir + os.sep + 'maise')

    # Get the Current Working Directory
    # cwd = os.getcwd()

    # Move to the actual directory where maise will run
    os.chdir(workdir)

    wf = open('maise.stdout', 'w')
    subprocess.call(['./maise'], stdout=wf)
    wf.close()

    if os.path.isfile('OSZICAR'):
        energies = np.loadtxt('OSZICAR')
    else:
        energies = None

    forces = None
    stress = None
    stress_kb = None
    if os.path.isfile('OUTCAR'):
        rf = open('OUTCAR', 'r')
        data = rf.read()

        pos_forces = re.findall(r'TOTAL-FORCE \(eV/Angst\)\s*-*\s*([-.\d\s]+)\s+-{2}', data)
        pos_forces = np.array([x.split() for x in pos_forces], dtype=float)

        if len(pos_forces) > 0 and len(pos_forces[-1]) % 7 == 0:
            pos_forces.shape = (len(pos_forces), -1, 7)
            forces = pos_forces[:, :, 3:6]
            # positions = pos_forces[:, :, :3]
        else:
            print('Forces and Positions could not be parsed : ', pos_forces.shape)
            print('pos_forces =\n%s ' % pos_forces)

        str_stress = re.findall('Total([.\d\s-]*)in', data)
        if len(str_stress) == 2:
            stress = np.array([[float(y) for y in x.split()] for x in str_stress])
        str_stress = re.findall('in kB([.\d\s-]*)energy', data)
        if len(str_stress) == 2:
            stress_kb = np.array([[float(y) for y in x.split()] for x in str_stress])

    create_new = False
    if not os.path.isfile('CONTCAR') or os.path.getsize("CONTCAR") == 0:
        create_new = True
        print('CONTCAR not found in entry: %s' % str(entry_id))
        i = 1
        while True:
            if not os.path.isfile('POSCAR-failed-%03s' % str(i)):
                os.rename('POSCAR', 'POSCAR-failed-%03s' % str(i))
                break
            else:
                i += 1
    else:
        new_structure = read_poscar('CONTCAR')
        # min_dist = np.min(new_structure.distance_matrix+np.ones((new_structure.natom,new_structure.natom)))
    min_dist = np.min(new_structure.distance_matrix()+(np.eye(new_structure.natom)*5))   # WIH
    print('Minimal distance= %8.7f' % min_dist)   # WIH

    if min_dist < 2.0:
        print('ERROR: MAISE finished with and structure with distances too close:', entry_id)  # WIH
        write_poscar(new_structure, workdir + os.sep + 'Collapsed_CONTCAR')  # WIH
        create_new = True   # WIH

    if create_new:
        new_structure = Structure.random_cell(structure.composition)
        ncalls = 0    # WIH

    if ncalls > max_ncalls:
        print('WARNING: Too many calls to MAISE and no relaxation succeeded, replacing structure: ', entry_id)    # WIH
        new_structure = Structure.random_cell(structure.composition)
        pcdb.entries.update({'_id': entry_id}, {'$set': {'status.ncalls': 0}})
        create_new = True
    else:
        pcdb.entries.update({'_id': entry_id}, {'$set': {'status.ncalls': ncalls}})
    pcdb.update(entry_id, structure=new_structure, properties={})

    # if not create_new and energies is not None and forces is not None and stress is not None:
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


def worker_vasp(db_settings, entry_id, workdir, relaxator_params):
    pcdb = get_database(db_settings)
    target_forces = relaxator_params['target_forces']
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


def worker_dftb(db_settings, entry_id, workdir, target_forces, relaxator_params):
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


def worker(db_settings, entry_id, workdir, target_forces, relaxator_params):
    pcdb = get_database(db_settings)
    pcm_log.info('[%s]: Starting relaxation. Target forces: %7.3e' % (str(entry_id), target_forces))

    if pcdb.is_locked(entry_id):
        return
    else:
        pcdb.lock(entry_id)
    structure = pcdb.get_structure(entry_id)
    structure = structure.scale()
    print('relaxator_params', relaxator_params)
    relaxer = IonRelaxation2(structure, workdir=workdir, target_forces=target_forces, waiting=False,
                             binary=relaxator_params['binary'], encut=1.3, kp_grid=None, kp_density=1E4,
                             relax_cell=True)
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


def is_evaluated(pcdb, entry_id, relaxator_params):
    status = get_current_status(pcdb, entry_id, relaxator_params)
    return status < relaxator_params['target_forces']


def get_current_status(pcdb, entry_id, relaxator_params, verbose=False):
    max_force = 1
    max_diag_stress = 1
    max_nondiag_stress = 1
    entry = pcdb.get_entry(entry_id)

    if entry is not None and 'properties' in entry and entry['properties'] is not None:
        if 'forces' in entry['properties'] and entry['properties']['forces'] is not None:
            forces = np.array(entry['properties']['forces']).reshape((-1, 3))
            max_force = np.max(np.apply_along_axis(np.linalg.norm, 1, forces))
        if 'stress' in entry['properties'] and entry['properties']['stress'] is not None:
            stress = np.array(entry['properties']['stress']).reshape((-1, 3))
            assert (stress.shape == (2, 3))
            max_diag_stress = np.abs(np.max(np.abs(stress[0])) - (relaxator_params['pressure'] / 1602.1766208))
            max_nondiag_stress = np.max(np.abs(stress[1]))
    else:
        pcm_log.debug('Bad entry')
        print(entry)
    if max(max_force, max_diag_stress, max_nondiag_stress) > relaxator_params['target_forces'] and verbose:
        if (max_force * max_diag_stress * max_nondiag_stress) == 1.0:
            print('No forces/stress information for entry: %s' % entry['_id'])
        else:
            print('Convergence status for entry: %s' % entry['_id'])
            print('Max Interatomic Force: %9.2E [eV/Ang]' % max_force)
            print('Max Diag stress      : %9.2E with pressure: %9.2E [eV/Ang^3]' % (max_diag_stress,
                                                                                    relaxator_params['pressure'] /
                                                                                    1602.1766208))
            print('Max NonDiag stress   : %9.2E [eV/Ang^3]' % max_nondiag_stress)
    return max(max_force, max_diag_stress, max_nondiag_stress)


if __name__ == '__main__':

    logging.basicConfig(level=logging.DEBUG)
    logger = logging.getLogger('pychemia')
    logger.addHandler(logging.NullHandler())
    logger.setLevel(logging.DEBUG)

    description = """Launch ab-initio codes for non-evaluated entries in a PyChemia Database"""

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
                        default=None, metavar='dbname', type=str, nargs='+',
                        help='PyChemia Database name (default: None)')
    parser.add_argument('-b', '--binary',
                        default='vasp', metavar='path', type=str,
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
    parser.add_argument('-s', '--source_dir',
                        default=None, metavar='path', type=str,
                        help='Working Directory  (default: None)')
    parser.add_argument('-l', '--slater_path',
                        default=None, metavar='path', type=str,
                        help='Slater Path, used only for DFTB+ (default: None)')
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
    db_settings = {'host': args.host, 'port': args.port, 'ssl': args.ssl, 'replicaset': args.replicaset}
    if args.user is not None:
        if args.passwd is None:
            raise ValueError('Password is mandatory if user is entered')
        db_settings['user'] = args.user
        db_settings['passwd'] = args.passwd

    relaxator_params = {'binary': args.binary,
                        'nmpiparal': args.nmpiparal,
                        'target_forces': args.target_forces,
                        'pressure': args.pressure_kB,
                        'source_dir': args.source_dir}

    print('pyChemia Evaluator using VASP')
    print('dbname    : %s' % args.dbname)
    print('source_dir: %s' % args.source_dir)
    print('host      : %s' % args.host)
    print('port      : %d' % args.port)
    print('user      : %s' % args.user)
    print('replicaset: %s' % args.replicaset)
    print('nparal    : %d' % args.nparal)
    print('nmpiparal : %d' % args.nmpiparal)
    print('binary    : %s' % str(args.binary))
    print('target-forces : %.2E' % args.target_forces)
    print('pressure_kB   : %.2E' % args.pressure_kB)
    print('evaluate_all  : %s' % str(args.evaluate_all))
    print('waiting       : %s' % str(args.waiting))
    print('ssl           : %s' % str(args.ssl))
    print('slater-path   : %s' % str(args.slater_path))

    print(db_settings)
    print(relaxator_params)

    worker = None
    if args.binary[:4].lower() == 'mais':
        worker = worker_maise
    elif args.binary[:4].lower() == 'dftb':
            worker = worker_dftb
    else:
        worker = worker_vasp

    evaluator = DirectEvaluator(db_settings=db_settings,
                                dbnames=args.dbname,
                                source_dir=args.source_dir,
                                is_evaluated=is_evaluated,
                                worker=worker,
                                worker_args=relaxator_params,
                                nconcurrent=args.nparal,
                                evaluate_all=args.evaluate_all)

    evaluator.run()
