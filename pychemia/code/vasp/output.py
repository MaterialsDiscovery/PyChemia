import re
import os
import numpy as np
from pychemia.utils.serializer import generic_serializer
from pychemia import pcm_log
from ..codes import CodeOutput
from .xml_output import parse_vasprun


class VaspOutput(CodeOutput):

    def __init__(self, filename='OUTCAR'):

        CodeOutput.__init__(self)
        if not os.path.isfile(filename):
            raise ValueError('File not found ' + filename)
        else:
            self.filename = filename
            rf = open(self.filename, 'r')
            self.data = rf.read()
            rf.close()

        self.charge = {}
        self.energies = None
        self.last_energy = None
        self.forces = None
        self.stress = None
        self.fermi = None
        self.bands = None
        self.positions = None
        self.kpoints = None
        self.array_sizes = {}
        self.species = None
        self.final_data = {}
        self.iteration_data = []
        self.outcar_parser()
        self.read()

    def read(self):
        if self.filename[-3:] == 'xml':
            ret = parse_vasprun(self.filename)
        else:
            self.read_outputfile(self.filename)
            self.outcar_parser()
            ret = {'energies': self.energies, 'forces': self.forces, 'bands': self.bands, 'positions': self.positions,
                   'kpoints': self.kpoints, 'stress': self.stress}
        self.output_values = ret

    def is_loaded(self):
        if self.output_values is None:
            return False
        else:
            return True

    def read_outputfile(self, filename=None):
        if filename is not None and os.path.exists(filename):
            rf = open(filename)
        else:
            rf = open(self.filename)
        self.data = rf.read()
        rf.close()

    def reload(self):
        rf = open(self.filename)
        self.data = rf.read()
        rf.close()

    def outcar_parser(self):

        for istr in ['NKPTS', 'NBANDS', 'NEDOS', 'NIONS', 'NGX', 'NGY', 'NGZ', 'NGXF', 'NGYF', 'NGZF', 'ISPIN']:
            redata = re.findall(istr + r'\s*=\s*(\d+)', self.data)
            if len(redata) > 0:
                self.array_sizes[istr] = int(redata[0])
        # pcm_log.info('Array sizes : ' + str(self.array_sizes))

        self.species = re.findall(r'POTCAR\s*:\s*[\w_]+\s*(\w+)', self.data)
        self.species = self.species[:int(len(self.species) / 2)]
        # pcm_log.info('Number of species (= number of POTCARs):' + str(self.species))

        pos_forces = re.findall(r'TOTAL-FORCE \(eV/Angst\)\s*-*\s*([-.\d\s]+)\s+-{2}', self.data)
        pos_forces = np.array([x.split() for x in pos_forces], dtype=float)

        if len(pos_forces) > 0 and len(pos_forces[-1]) % 6 == 0:
            pos_forces.shape = (len(pos_forces), -1, 6)
            forces = pos_forces[:, :, 3:]
            positions = pos_forces[:, :, :3]
            pcm_log.debug('Positions from OUTCAR: %d iterations' % len(positions))
            pcm_log.debug('Forces from OUTCAR: %d iterations' % len(forces))

            self.forces = forces
            self.positions = positions
            self.array_sizes['NIONSTEPS'] = len(self.forces)
            pcm_log.debug('Number of Ionic steps: ' + str(self.array_sizes['NIONSTEPS']))
        else:
            print('Forces and Positions could not be parsed : ', pos_forces.shape)
            print('pos_forces =\n%s ' % pos_forces)

        fermi = re.findall(r'E-fermi\s+:\s+([-.\d]+)', self.data)
        fermi = np.array(fermi, dtype=float)
        # pcm_log.debug('Fermi Level(eV): ' + str(fermi))
        self.fermi = fermi

        # This regex covers the entire information for each electronic iteration.
        # The expression in the middle, catch everything but ('>' greater than)
        # The final part catch the value of the energy
        # It returns a list of tuples in the form [(ionic iter, Energy (elect. Iter), ...]
        energy = re.findall(r'Iteration\s*(\d+)\s*\(\s*(\d+)\)[-+*/=():.\s\d\w]+>0\)\s*=\s*([-.\d]+)', self.data)
        # pcm_log.debug('Energy(eV) [(ionic iter, Energy (elect. Iter)),...]: ' + str(energy))

        self.energies = []
        index = None

        for ienergy in energy:
            self.energies.append([int(ienergy[0]), int(ienergy[1]), float(ienergy[2])])

        ### BAD idea this way
        ###re.findall(r'free\s*energy\s*TOTEN\s*=\s*([-.\d]+)\s*eV', data)
        if not energy:
            energy = re.findall(
                r'FREE ENERGIE OF THE ION-ELECTRON SYSTEM\s*\(eV\)\s*[-]+\s*free\s*energy\s*TOTEN\s*=\s*([-.\d]+)\s*eV',
                self.data)

            counter = 0
            for ienergy in energy:
                self.energies.append([counter, 1, float(ienergy)])
                counter += 1

        level0 = 0
        level1 = 0
        for i in self.energies:
            if i[0] > level0:
                self.last_energy = i[2]
                level0 = i[0]
            elif i[0] == level0 and i[1] > level1:
                self.last_energy = i[2]
                level1 = i[1]

        pattern = r"k-point([\s\d]+):([\d\s.+-]+)\s*band No.\s*band energies\s*occupation\s*\n([\d\w\s.+-]+)\n\n"
        bands = re.findall(pattern, self.data)
        bands_dict = {}
        for iband in bands:
            try:
                ikpt = int(iband[0])
                kpt_pos = [float(x) for x in iband[1].split()]
                kpt_eigenvals = np.array([float(x) for x in iband[2].split()[:3 * self.array_sizes['NBANDS']]])
                kpt_eigenvals.reshape((-1, 3))
                bands_dict[ikpt] = {'position': kpt_pos, 'eigenvalues': kpt_eigenvals}
            except ValueError:
                print('Error parsing bands')
        self.bands = bands_dict

        if False and 'NIONSTEPS' in self.array_sizes:
            if len(bands) != self.array_sizes['NBANDS'] * self.array_sizes['ISPIN'] * self.array_sizes['NKPTS'] \
                    * self.array_sizes['NIONSTEPS']:
                pcm_log.debug('NBANDS: %s != ISPIN: %s x NKPTS: %s x NIONSTEPS: %s' % (self.array_sizes['NBANDS'],
                                                                                       self.array_sizes['ISPIN'],
                                                                                       self.array_sizes['NKPTS'],
                                                                                       self.array_sizes['NIONSTEPS']))
        # pcm_log.info('Bands : ' + str(bands))

        stress = re.findall(r'in\s+kB ([-*.\s\d]+)external', self.data)
        # Converted from kBar to GPa
        if stress:
            stress = 0.1 * np.array([x.split() if '*' not in x else 6 * [float('nan')] for x in stress], dtype=float)
            stress.reshape((-1, 6))
            nionsteps = len(stress)
            self.stress = np.zeros((nionsteps, 3, 3))
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

        charge = re.findall(r'total\s*charge\s*#\s*of\s*ion\s*s\s*p\s*d\s*tot\s*-+([-.\s\d]+)\s--', self.data)
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

        # Processing Final
        self.free_energy()
        self.iterations()

    def free_energy(self):
        subdata = re.findall('FREE ENERGIE OF THE ION-ELECTRON SYSTEM [\w\d\s()\n-=>]*\n\s*\n\n\n-', self.data)
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

    @property
    def energy(self):
        if len(self.energies) > 0:
            return self.energies[-1][-1]

    def iterations(self):
        # Capture the data for each iteration. The symbol '>' is not included on purpose to close
        # the capture close to "energy(sigma->0)..."
        subdata = re.findall(r'-+ Iteration\s*(\d+)\(\s*(\d+)\)\s*-+\s*([\s\d\w:.,\-=()/+*]*)energy', self.data)
        for i in subdata:
            io_iter = int(i[0]) - 1
            scf_iter = int(i[1]) - 1
            # print io_iter, scf_iter
            while len(self.iteration_data) <= io_iter:
                self.iteration_data.append([])
            while len(self.iteration_data[io_iter]) <= scf_iter:
                self.iteration_data[io_iter].append([])
            self.iteration_data[io_iter][scf_iter] = {}
            self.iteration_data[io_iter][scf_iter]['Timing'] = {}
            for j in ['POTLOK', 'SETDIJ', 'EDDIAG', 'RMM-DIIS', 'ORTHCH', 'DOS', 'CHARGE', 'MIXING', 'LOOP']:
                iter_block = re.findall(j + r':([\d\w:. ]*)\n', i[2])
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
                iter_block = re.findall(j + r'\s*=\s*([\d.-]*)[\s\w\d]*\n', i[2])
                if len(iter_block) > 0:
                    self.iteration_data[io_iter][scf_iter]['Free_Energy'][j] = float(iter_block[0])

    def has_forces_stress_energy(self):
        return self.forces is not None and self.stress is not None and self.last_energy is not None

    def __str__(self):
        ret = '\nForces:\n'
        index = 0
        for iforce in self.forces[-1]:
            index += 1
            ret += "%3d %12.6f %12.6f %12.6f\n" % (index, iforce[0], iforce[1], iforce[2])
        ret += '\nStress:\n'
        for j in range(3):
            ret += '    %12.6f %12.6f %12.6f\n' % tuple(self.stress[-1][j])
        ret += '\n'

        ret += 'Free Energy: %12.6f\n' % self.energy
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
        for i in ['charge', 'energy', 'forces', 'stress']:
            ret[i] = eval('self.' + i)

        for i in ret:
            if isinstance(ret[i], np.ndarray):
                ret[i] = generic_serializer(ret[i])
        return ret

    def get_magnetization(self):
        ret = {}
        for mag_dir in ['x', 'y', 'z']:
            mag = re.findall(r'magnetization\s*\(%s\)\s*#\s*of\s*ion\s*s\s*p\s*d\s*tot\s*-+([-.\s\d]+)\s--' % mag_dir,
                             self.data)
            if mag:
                mag = np.array([x.split() for x in mag], dtype=float)
                mag.shape = (len(mag), -1, 5)
                if len(mag) == 2:
                    assert np.all(np.array(mag[0] == mag[1]))
                mag = mag[-1, :, 1:]
                ret[mag_dir] = {}
                ret[mag_dir]['s'] = [float(x) for x in mag[:, 0]]
                ret[mag_dir]['p'] = [float(x) for x in mag[:, 1]]
                ret[mag_dir]['d'] = [float(x) for x in mag[:, 2]]
                ret[mag_dir]['total'] = [float(x) for x in mag[:, 3]]
        return ret

    def get_memory_used(self):

        ret = {}
        datablock = re.findall(r"total amount of memory used by VASP on root node([\w\s\d.=:-]*)\n \n {2}", self.data)
        if len(datablock) != 1:
            return None
        else:
            datablock = datablock[0]
        print(datablock)

        for iline in datablock.split('\n'):
            if ':' in iline:
                key_value = iline.split(':')
                try:
                    ret[key_value[0].strip()] = (float(key_value[1].split()[0]), key_value[1].split()[1])
                except ValueError:
                    print('Failed to parse: ', key_value)
        return ret

    def get_general_timing(self):
        ret = {}
        last_lines = self.data.split('\n')[-15:]
        if 'General timing and accounting informations for this job' in last_lines[0]:
            for iline in last_lines[2:]:
                if ':' in iline:
                    key_value = iline.split(':')
                    try:
                        ret[key_value[0].strip()] = float(key_value[1].strip())
                    except ValueError:
                        print('Failed to parse: ', key_value)
        return ret

    @property
    def is_finished(self):
        return len(self.get_general_timing()) > 0
