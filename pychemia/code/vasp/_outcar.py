__author__ = 'Guillermo Avendano-Franco'

import os


class VaspOutput():

    def __init__(self, filename='OUTCAR'):

        self.magnetization = None
        self.total_charge = None
        self.free_energy = None
        self.forces = None
        self.stress = None

        if not os.path.isfile(filename):
            raise ValueError('File not found '+filename)

        rf = open(filename)
        lines = rf.readlines()

        for i in range(len(lines)):
            if lines[i].strip() == 'magnetization (x)':
                self.get_magnetization(i, lines)
            elif lines[i].strip() == 'total charge':
                self.get_total_charge(i, lines)
            elif lines[i].strip() == 'FREE ENERGIE OF THE ION-ELECTRON SYSTEM (eV)':
                self.get_free_energy(i, lines)
            elif lines[i].strip() == 'FORCE on cell =-STRESS in cart. coord.  units (eV):':
                self.get_stress(i, lines)
            elif lines[i].strip() == 'POSITION                                       TOTAL-FORCE (eV/Angst)':
                self.get_forces(i, lines)

    def get_magnetization(self, iline, lines):
        keys = [x.strip() for x in lines[iline+2].split(5*' ')]
        i = 3
        self.magnetization = {}
        for k in keys:
            self.magnetization[k] = []
        while True:
            i += 1
            if lines[iline+i].strip().startswith('-'):
                break
            values = [float(x) for x in lines[iline + i].split()]
            assert(len(keys) == len(values))
            for ik in range(len(keys)):
                self.magnetization[keys[ik]].append(values[ik])

    def get_total_charge(self, iline, lines):
        keys = [x.strip() for x in lines[iline+2].split(5*' ')]
        i = 3
        self.total_charge = {}
        for k in keys:
            self.total_charge[k] = []
        while True:
            i += 1
            if lines[iline+i].strip().startswith('-'):
                break
            values = lines[iline + i].split()
            assert(len(keys) == len(values))
            for ik in range(len(keys)):
                self.total_charge[keys[ik]].append(values[ik])

    def get_free_energy(self, iline, lines):
        self.free_energy = float(lines[iline+2].split()[4])

    def get_forces(self, iline, lines):
        i = 1
        self.forces = []
        while True:
            i += 1
            if lines[iline+i].strip().startswith('-'):
                break
            values = lines[iline + i].split()
            self.forces.append([float(x) for x in values[-3:]])

    def get_stress(self, iline, lines):
        i = 2
        self.stress = {}

        while True:
            i += 1
            if lines[iline+i].strip().startswith('-'):
                break
            key = lines[iline + i][:9].strip()
            #print lines[iline + i][9:]
            values = [float(x) for x in lines[iline + i][9:].split()]
            #print values
            self.stress[key] = values
        i += 1
        values = [float(x) for x in lines[iline + i][9:].split()]
        self.stress['Total'] = values

    def __str__(self):
        ret = '\nForces:\n'
        index = 0
        for iforce in self.forces:
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
