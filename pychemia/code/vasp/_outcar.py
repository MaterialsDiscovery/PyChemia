__author__ = 'Guillermo Avendano-Franco'

import re
import os
import logging
import numpy as np

logging.basicConfig(level=logging.INFO)


class VaspOutput():
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

        if not os.path.isfile(filename):
            raise ValueError('File not found ' + filename)

    def outcar_parser(self):
        rf = open(self.filename, 'r')
        data = rf.read()
        rf.close()

        for istr in ['NKPTS', 'NBANDS', 'NEDOS', 'NIONS', 'NGX', 'NGY', 'NGZ', 'NGXF', 'NGYF', 'NGZF', 'ISPIN']:
            self.array_sizes[istr] = int(re.findall(istr + r'\s*=\s*(\d+)', data)[0])
        logging.info('Array sizes : ' + str(self.array_sizes))

        self.species = re.findall(r'POTCAR\s*:\s*[\w_]+\s*(\w+)', data)
        self.species = self.species[:len(self.species) / 2]
        logging.info('Number of species (= number of POTCARs):' + str(self.species))

        pos_forces = re.findall(r'TOTAL-FORCE \(eV/Angst\)\s*-*\s*([-.\d\s]+)\s+-{2}', data)
        pos_forces = np.array([x.split() for x in pos_forces], dtype=float)
        pos_forces.shape = (len(pos_forces), -1, 6)
        forces = pos_forces[:, :, 3:]
        positions = pos_forces[:, :, :3]
        logging.debug('Positions from OUTCAR: \n' + str(positions))
        logging.debug('Forces from OUTCAR: \n' + str(forces))

        self.forces = forces
        self.positions = positions
        self.array_sizes['NIONSTEPS'] = len(self.forces)
        logging.info('Number of Ionic steps: ' + str(self.array_sizes['NIONSTEPS']))

        fermi = re.findall(r'E-fermi\s+:\s+([-.\d]+)', data)
        fermi = np.array(fermi, dtype=float)
        logging.debug('Fermi Level(eV): ' + str(fermi))
        self.fermi = fermi

        # This regex covers the entire information for each electronic iteration.
        # The expression in the middle, catch everything but ('>' greater than)
        # The final part catch the value of the energy
        # It returns a list of tuples in the form [(ionic iter, Energy (elect. Iter), ...]
        energy = re.findall(r'Iteration\s*(\d+)[-+*/=():.\s\d\w]+>0\)\s*=\s*([-.\d]+)', data)
        logging.debug('Energy(eV) [(ionic iter, Energy (elect. Iter)),...]: ' + str(energy))
        self.energy = energy

        # TODO: Check for magnetic cases
        kpoints = re.findall(r'k-point\s*\d+\s*:\s*([-.\d\s]+)band', data)
        kpoints = np.array([x.split() for x in kpoints], dtype=float)

        if len(kpoints) == 2 * self.array_sizes['NKPTS'] * self.array_sizes['NIONSTEPS']:
            assert (self.array_sizes['ISPIN'] == 2)
        else:
            assert (self.array_sizes['ISPIN'] == 1)

        kpoints.shape = (self.array_sizes['NIONSTEPS'], self.array_sizes['ISPIN'], self.array_sizes['NKPTS'], 3)
        logging.debug('K-points:' + str(kpoints))
        logging.debug('K-points shape;:' + str(kpoints.shape))

        if self.array_sizes['ISPIN'] == 2:
            assert (np.all(kpoints[:, 0, :, :] == kpoints[:, 1, :, :]))
        kpoints = kpoints[:, 0, :, :]
        self.kpoints = kpoints

        bands = re.findall(r'^\s*[1-9]\d*\s+([-\d]+\.[\d]+)\s+[-.\d]+\s*\n', data, re.MULTILINE)
        assert (len(bands) == self.array_sizes['NBANDS'] * self.array_sizes['ISPIN']
                * self.array_sizes['NKPTS']
                * self.array_sizes['NIONSTEPS'])
        bands = np.array(bands, dtype=float)
        logging.info('Bands : ' + str(bands))

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
            logging.info('Stress (GPa):\n ' + str(self.stress))

        charge = re.findall(r'total\s*charge\s*#\s*of\s*ion\s*s\s*p\s*d\s*tot\s*-+([-.\s\d]+)\s--', data)
        if charge:
            charge = np.array([x.split() for x in charge], dtype=float)
            charge.shape = (len(charge), -1, 5)
            if len(charge) == 2:
                assert np.all(charge[0] == charge[1])
            charge = charge[0, :, 1:]
            self.charge['s'] = charge[:, 0]
            self.charge['p'] = charge[:, 1]
            self.charge['d'] = charge[:, 2]
            self.charge['total'] = charge[:, 3]
            logging.debug('Charge :' + str(self.charge))

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
            logging.debug('Magnetization :' + str(self.magnetization))

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

    def to_dict(self):
        ret = {}
        for i in ['magnetization', 'charge', 'energy', 'forces', 'stress']:
            ret[i] = eval('self.' + i)
        return ret
