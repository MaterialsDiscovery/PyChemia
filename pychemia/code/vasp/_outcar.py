__author__ = 'Guillermo Avendano-Franco'

import re
import os
import numpy as np
from pychemia.serializer import generic_serializer
from pychemia import pcm_log


class VaspOutput:
    def __init__(self, filename='OUTCAR'):

        self.magnetization = {}
        self.charge = {}
        self.energy = None
        self.forces = None
        self.stress = None
        self.fermi = None
        self.positions = None
        self.kpoints = None
        self.filename = filename
        self.array_sizes = {}
        self.species = None
        self.final_data = {}
        self.iteration_data = []

        if not os.path.isfile(filename):
            raise ValueError('File not found ' + filename)
        else:
            self.outcar_parser()

    def outcar_parser(self):

        rf = open(self.filename, 'r')
        data = rf.read()
        rf.close()

        for istr in ['NKPTS', 'NBANDS', 'NEDOS', 'NIONS', 'NGX', 'NGY', 'NGZ', 'NGXF', 'NGYF', 'NGZF', 'ISPIN']:
            redata = re.findall(istr + r'\s*=\s*(\d+)', data)
            if len(redata) > 0:
                self.array_sizes[istr] = int(redata[0])
        # pcm_log.info('Array sizes : ' + str(self.array_sizes))

        self.species = re.findall(r'POTCAR\s*:\s*[\w_]+\s*(\w+)', data)
        self.species = self.species[:len(self.species) / 2]
        # pcm_log.info('Number of species (= number of POTCARs):' + str(self.species))

        pos_forces = re.findall(r'TOTAL-FORCE \(eV/Angst\)\s*-*\s*([-.\d\s]+)\s+-{2}', data)
        pos_forces = np.array([x.split() for x in pos_forces], dtype=float)

        if len(pos_forces) > 0 and len(pos_forces[-1]) % 6 == 0:
            pos_forces.shape = (len(pos_forces), -1, 6)
            forces = pos_forces[:, :, 3:]
            positions = pos_forces[:, :, :3]
            # pcm_log.debug('Positions from OUTCAR: \n' + str(positions))
            # pcm_log.debug('Forces from OUTCAR: \n' + str(forces))

            self.forces = forces
            self.positions = positions
            self.array_sizes['NIONSTEPS'] = len(self.forces)
            # pcm_log.info('Number of Ionic steps: ' + str(self.array_sizes['NIONSTEPS']))
        else:
            print 'Forces and Positions could not be parsed : ', pos_forces.shape

        fermi = re.findall(r'E-fermi\s+:\s+([-.\d]+)', data)
        fermi = np.array(fermi, dtype=float)
        # pcm_log.debug('Fermi Level(eV): ' + str(fermi))
        self.fermi = fermi

        # This regex covers the entire information for each electronic iteration.
        # The expression in the middle, catch everything but ('>' greater than)
        # The final part catch the value of the energy
        # It returns a list of tuples in the form [(ionic iter, Energy (elect. Iter), ...]
        energy = re.findall(r'Iteration\s*(\d+)[-+*/=():.\s\d\w]+>0\)\s*=\s*([-.\d]+)', data)
        # pcm_log.debug('Energy(eV) [(ionic iter, Energy (elect. Iter)),...]: ' + str(energy))

        self.energy = []
        index = None

        for i in energy:
            if index != i[0]:
                if len(self.energy) == int(i[0]) - 1:
                        self.energy.append([])
                self.energy[int(i[0])-1].append(float(i[1]))

        # TODO: Check for magnetic cases
        kpoints = re.findall(r'k-point\s*\d+\s*:\s*([-.\d\s]+)band', data)
        kpoints = np.array([x.split() for x in kpoints], dtype=float)

        if False and self.array_sizes['NKPTS'] < 1000:
            if len(kpoints) == 2 * self.array_sizes['NKPTS'] * self.array_sizes['NIONSTEPS']:
                assert (self.array_sizes['ISPIN'] == 2)
            else:
                assert (self.array_sizes['ISPIN'] == 1)

            kpoints.shape = (self.array_sizes['NIONSTEPS'], self.array_sizes['ISPIN'], self.array_sizes['NKPTS'], 3)
            # pcm_log.debug('K-points:' + str(kpoints))
            # pcm_log.debug('K-points shape;:' + str(kpoints.shape))

            if self.array_sizes['ISPIN'] == 2:
                assert (np.all(kpoints[:, 0, :, :] == kpoints[:, 1, :, :]))
            kpoints = kpoints[:, 0, :, :]
            self.kpoints = kpoints

        bands = re.findall(r'^\s*[1-9]\d*\s+([-\d]+\.[\d]+)\s+[-.\d]+\s*\n', data, re.MULTILINE)

        if 'NIONSTEPS' in self.array_sizes:
            if len(bands) != self.array_sizes['NBANDS'] * self.array_sizes['ISPIN'] * self.array_sizes['NKPTS'] \
                    * self.array_sizes['NIONSTEPS']:
                pcm_log.error('Bad number of bands')
        bands = np.array(bands, dtype=float)
        # pcm_log.info('Bands : ' + str(bands))

        stress = re.findall(r'in\s+kB ([-.\s\d]+)external', data)
        # Converted from kBar to GPa
        if stress:
            stress = 0.1 * np.array([x.split() for x in stress], dtype=float)
            stress.shape = (self.array_sizes['NIONSTEPS'], 6)
            self.stress = np.zeros((self.array_sizes['NIONSTEPS'], 3, 3))
            self.stress[:, 0, 0] = stress[:, 0]
            self.stress[:, 1, 1] = stress[:, 1]
            self.stress[:, 2, 2] = stress[:, 2]
            self.stress[:, 0, 1] = stress[:, 3]
            self.stress[:, 1, 2] = stress[:, 4]
            self.stress[:, 0, 2] = stress[:, 5]
            self.stress[:, 1, 0] = stress[:, 3]
            self.stress[:, 2, 1] = stress[:, 4]
            self.stress[:, 2, 0] = stress[:, 5]
            # pcm_log.info('Stress (GPa):\n ' + str(self.stress))

        charge = re.findall(r'total\s*charge\s*#\s*of\s*ion\s*s\s*p\s*d\s*tot\s*-+([-.\s\d]+)\s--', data)
        if charge:
            charge = np.array([x.split() for x in charge], dtype=float)
            charge.shape = (len(charge), -1, 5)
            if len(charge) == 2:
                if np.all(charge[0] != charge[1]):
                    pcm_log.error('Bad Charge')
            charge = charge[0, :, 1:]
            self.charge['s'] = list(charge[:, 0])
            self.charge['p'] = list(charge[:, 1])
            self.charge['d'] = list(charge[:, 2])
            self.charge['total'] = list(charge[:, 3])
            # pcm_log.debug('Charge :' + str(self.charge))

        magnetization = re.findall(r'magnetization\s*\(x\)\s*#\s*of\s*ion\s*s\s*p\s*d\s*tot\s*-+([-.\s\d]+)\s--', data)
        if magnetization:
            magnetization = np.array([x.split() for x in magnetization], dtype=float)
            magnetization.shape = (len(magnetization), -1, 5)
            if len(magnetization) == 2:
                assert np.all(magnetization[0] == magnetization[1])
            magnetization = magnetization[0, :, 1:]
            self.magnetization['s'] = magnetization[:, 0]
            self.magnetization['p'] = magnetization[:, 1]
            self.magnetization['d'] = magnetization[:, 2]
            self.magnetization['total'] = magnetization[:, 3]
            # pcm_log.debug('Magnetization :' + str(self.magnetization))

        # Processing Final Data
        self.free_energy(data)
        self.iterations(data)

    def free_energy(self, data):
        subdata = re.findall('FREE ENERGIE OF THE ION-ELECTRON SYSTEM [\w\d\s\(\)\n-=>]*\n \n\n\n-', data)
        if len(subdata) == 1:
            data_split = subdata[0].split()
            try:
                float(data_split[12])
                float(data_split[17])
                float(data_split[20])
            except ValueError:
                raise ValueError('Energy values could not be parsed from: %s' % subdata)
            self.final_data['energy'] = {'free_energy': float(data_split[12]),
                                         'energy_without_entropy': float(data_split[17]),
                                         'energy(sigma->0)': float(data_split[20])}

    def iterations(self, data):
        # Capture the data for each iteration. The symbol '>' is not included on purpose to close
        # the capture close to "energy(sigma->0)..."
        subdata = re.findall(r'-+ Iteration\s*(\d+)\(\s*(\d+)\)\s*-+\s*([\s\d\w:.,\-=()/+\*]*)energy', data)
        for i in subdata:
            io_iter = int(i[0])-1
            scf_iter = int(i[1])-1
            # print io_iter, scf_iter
            while len(self.iteration_data) <= io_iter:
                self.iteration_data.append([])
            while len(self.iteration_data[io_iter]) <= scf_iter:
                self.iteration_data[io_iter].append([])
            self.iteration_data[io_iter][scf_iter] = {}
            self.iteration_data[io_iter][scf_iter]['Timing'] = {}
            for j in ['POTLOK', 'SETDIJ', 'EDDIAG', 'RMM-DIIS', 'ORTHCH', 'DOS', 'CHARGE', 'MIXING', 'LOOP']:
                iter_block = re.findall(j+r':([\d\w:. ]*)\n', i[2])
                if len(iter_block) > 0:
                    iter_block_split = iter_block[0].split()
                    try:
                        float(iter_block_split[2][:-1])
                        float(iter_block_split[5])
                    except ValueError:
                        raise ValueError('Could not retrieve data', iter_block)
                    self.iteration_data[io_iter][scf_iter]['Timing'][j] = {'cpu': float(iter_block_split[2][:-1]),
                                                                           'real': float(iter_block_split[5])}

            self.iteration_data[io_iter][scf_iter]['Free_Energy'] = {}
            for j in ['PSCENC', 'TEWEN', 'DENC', 'EXHF', 'XCENC', 'EENTRO', 'EBANDS', 'EATOM', 'TOTEN']:
                iter_block=re.findall(j+r'\s*=\s*([\d\.-]*)[\s\w\d]*\n', i[2])
                if len(iter_block) > 0:
                    self.iteration_data[io_iter][scf_iter]['Free_Energy'][j] = float(iter_block[0])

    def has_forces_stress_energy(self):
        return self.forces is not None and self.stress is not None and self.energy is not None

    def __str__(self):
        ret = '\nForces:\n'
        index = 0
        for iforce in self.forces[-1]:
            index += 1
            ret += "%3d %9.6f %9.6f %9.6f\n" % (index, iforce[0], iforce[1], iforce[2])
        ret += '\nStress:\n'
        for istress in sorted(self.stress):
            if istress != 'Total':
                ret += '%8s ' % istress
                for i in self.stress[istress]:
                    ret += ' %12.6f' % i
                ret += '\n'
        ret += '%8s ' % 'Total'
        for i in self.stress['Total']:
            ret += ' %12.6f' % i
        ret += '\n'

        ret += '\nFree Energy:\n'
        ret += str(self.free_energy)
        return ret

    def relaxation_info(self):

        info = {}
        if self.stress is not None:
            info['avg_stress_diag'] = np.average(np.abs(self.stress[-1].diagonal()))
            info['avg_stress_non_diag'] = np.average(np.abs(np.triu(self.stress[-1], 1)))
        if self.forces is not None:
            info['avg_force'] = np.average(np.abs(np.apply_along_axis(np.linalg.norm, 1, self.forces[-1])))
        return info

    @property
    def to_dict(self):
        ret = {}
        for i in ['magnetization', 'charge', 'energy', 'forces', 'stress']:
            ret[i] = eval('self.' + i)

        for i in ret:
            if isinstance(ret[i], np.ndarray):
                ret[i] = generic_serializer(ret[i])
        return ret


def read_vasp_stdout(filename):
    rf = open(filename)
    data = rf.read()

    re_str = r'\n([\w]{3})\:\s*([\d])+\s*([\dE+-.]+)\s*([\dE+-.]+)\s*([\dE+-.]+)\s*([\d]+)\s*([\dE+-.]+)\s*[-\ .\w+]+'

    result = re.findall(re_str, data)

    ret = []
    for i in result:
        ret.append([i[0], int(i[1]), float(i[2]), float(i[3]), float(i[4]), int(i[5]), float(i[6])])

    number_of_scf_per_ionic_iter = []
    final_energy_after_scf = []
    index = 0
    counter = 0
    for i in ret:
        if i[1] > index:
            index = i[1]
        else:
            number_of_scf_per_ionic_iter.append(index)
            final_energy_after_scf.append(ret[counter][2])
            index = 0
        counter += 1

    return {'iterations': number_of_scf_per_ionic_iter, 'energies': final_energy_after_scf, 'data': ret}
