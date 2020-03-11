import os
import re
import subprocess

import numpy as np

from pychemia import pcm_log, Structure
from pychemia.utils.serializer import generic_serializer
from pychemia.utils.periodic import atomic_number, atomic_symbol
from ..codes import CodeRun


class FireBall(CodeRun):
    def __init__(self, workdir='.', fdata_path=None):
        CodeRun.__init__(self, executable='fireball.x', workdir=workdir, use_mpi=False)
        self.workdir = None
        # The five sections of Fireball
        self.option = {}
        self.tds = {}
        self.quench = {}
        self.tds = {}
        self.mesh = {}
        self.output = {}

        self.fdata_path = os.path.abspath(fdata_path)
        self.structure = None
        self.runner = None
        self.kpoints = None

    @property
    def to_dict(self):
        ret = {}
        if self.option:
            ret['option'] = self.option
        if self.output:
            ret['output'] = self.output
        if self.quench:
            ret['quench'] = self.quench
        if self.mesh:
            ret['mesh'] = self.mesh
        if self.tds:
            ret['tds'] = self.tds
        return ret

    @staticmethod
    def from_dict(dictionary):
        ret = FireBall()
        if 'option' in dictionary:
            ret.option = dictionary['option']
        if 'output' in dictionary:
            ret.option = dictionary['output']
        if 'quench' in dictionary:
            ret.option = dictionary['quench']
        if 'mesh' in dictionary:
            ret.option = dictionary['mesh']
        if 'tds' in dictionary:
            ret.option = dictionary['tds']
        return ret

    def initialize(self, structure, workdir=None, kpoints=None, executable='fireball.x'):
        assert structure.is_perfect
        self.structure = structure
        if workdir is not None:
            self.workdir = workdir
        else:
            self.workdir = '.'
        if not os.path.lexists(self.workdir):
            os.mkdir(self.workdir)
        self.kpoints = kpoints
        self.executable = executable

    def set_inputs(self, rms=0.1):
        self.write_input(filename=self.workdir + os.sep + 'fireball.in')
        self.write_basis(filename=self.workdir + os.sep + 'input.bas')
        if self.structure.is_periodic:
            self.write_lattice(filename=self.workdir + os.sep + 'input.lvs')
            self.write_kpoints(filename=self.workdir + os.sep + 'input.kpts')
        self.link_fdata()
        wf = open(self.workdir + os.sep + 'rms.input', 'w')
        wf.write('%f\n' % rms)
        wf.close()

    def get_outputs(self):
        pass

    def run(self, num_threads=None, mpi_num_procs=None, nodefile=None, wait=True, verbose=False):
        cwd = os.getcwd()
        os.chdir(self.workdir)
        stdout = open('fireball.log', 'w')
        sp = subprocess.Popen(self.executable, stdout=stdout)
        os.chdir(cwd)
        self.runner = sp
        return sp

    def run_status(self):
        if self.runner is None:
            pcm_log.info('Fireball not finish')
            filename = self.workdir + os.sep + 'fireball.log'
            if os.path.exists(filename):
                read_fireball_stdout(filename=filename)
            return
        if self.runner.poll() == 0:
            pcm_log.info('Fireball complete normally')

    def finalize(self):
        pass

    def write_input(self, filename='fireball.in'):

        wf = open(filename, 'w')
        wf.write(str(self))
        wf.close()

    def __str__(self):
        ret = ''
        if self.option:
            ret += '&OPTION\n'
            for variable in sorted(self.option):
                if isinstance(self.option[variable], str):
                    ret += "%s = '%s'\n" % (variable, str(self.option[variable]))
                else:
                    ret += '%s = %s\n' % (variable, str(self.option[variable]))
            ret += '&END\n'

        if self.output:
            ret += '&OUTPUT\n'
            for variable in sorted(self.output):
                ret += '%s = %s\n' % (variable, str(self.output[variable]))
            ret += '&END\n'

        if self.quench:
            ret += '&QUENCH\n'
            for variable in sorted(self.quench):
                ret += '%s = %s\n' % (variable, str(self.quench[variable]))
            ret += '&END\n'

        if self.mesh:
            ret += '&MESH\n'
            for variable in sorted(self.mesh):
                ret += '%s = %s\n' % (variable, str(self.mesh[variable]))
            ret += '&END\n'

        if self.tds:
            ret += '&TDS\n'
            for variable in sorted(self.tds):
                ret += '%s = %s\n' % (variable, str(self.tds[variable]))
            ret += '&END\n'

        return ret

    def write_basis(self, filename='input.bas'):
        write_geometry_bas(self.structure, filename)

    def write_lattice(self, filename='input.lvs'):
        wf = open(filename, 'w')
        for i in range(3):
            wf.write('%13.7f %13.7f %13.7f\n' % (self.structure.cell[i, 0],
                                                 self.structure.cell[i, 1],
                                                 self.structure.cell[i, 2]))
        wf.close()

    def write_kpoints(self, filename):
        pass

    def link_fdata(self):
        if self.fdata_path is not None:
            linkpath = self.workdir + os.sep + 'Fdata'
            if os.path.exists(linkpath):
                os.remove(linkpath)
            os.symlink(self.fdata_path, linkpath)

    def cluster_relaxation(self):
        self.option = {'basisfile': 'input.bas', 'nstepf': 5000, 'iquench': -3, 'icluster': 1, 'iqout': 3, 'dt': 0.25}
        self.output = {'iwrtxyz': 1}


def read_fireball_stdout(filename):
    rf = open(filename, 'r')
    data = rf.read()

    atom_data = re.findall(r'Atom Coordinates from Basis File:([\d\w\s=\-.#]+)===\n', data)[0]

    symbols = []
    positions = []
    for iline in atom_data.split('\n'):
        # print iline
        tmp = iline.split()
        if len(tmp) == 6 and tmp[0].isdigit():
            symbols.append(tmp[1])
            positions.append([float(x) for x in tmp[2:5]])

    natom = len(symbols)
    initial_positions = np.array(positions)
    # print symbols
    # print initial_positions

    forces_data = re.findall(r'The grand total force \(eV/A\):([\w\d\s.\-=+]*)Cartesian', data)
    # print len(forces_data)

    forces = np.zeros((len(forces_data), natom, 3))
    # print forces.shape

    iteration = 0
    for idata in forces_data:
        for iline in idata.split('\n'):
            if 'iatom' in iline:
                fields = iline.split()
                forces[iteration, int(fields[2]) - 1] = np.array(fields[-3:], dtype=float)
        iteration += 1

    ret = re.findall('Cartesian Forces:\s*Max = ([=\w\s\d.]*)RMS = ([\s\d.]*)\n', data)
    ret = [[float(x) for x in y] for y in ret]
    max_force = [x[0] for x in ret]
    rms = [x[1] for x in ret]

    energy_data = re.findall(r'---------- T H E\s+T O T A L\s+E N E R G Y -----------([\s\w\d.\-=/]+)--- \n', data)

    energy = []
    for idata in energy_data:
        for iline in idata.split('\n'):
            if 'Time step' in iline:
                tmp = iline.split()
                ienergy = {'Time_step': int(tmp[3]), 'SCF_step': int(tmp[7]), 'etot/atom': float(tmp[10])}
            elif len(iline.split('=')) == 2:
                tmp = iline.split('=')
                ienergy[tmp[0].strip()] = float(tmp[1])
        energy.append(ienergy)

    gt_data = re.findall(r'Grand Total = nuclear kinetic \+ potential = [\s\-\d.]*', data)
    gt = [float(i.split()[-1]) for i in gt_data]

    gtpa_data = re.findall(r'grand total energy per atom = [\s\-\d.]*', data)
    gtpa = [float(i.split()[-1]) for i in gtpa_data]

    ret = {'symbols': symbols,
           'initial_positions': initial_positions,
           'forces': generic_serializer(forces),
           'energetics': energy,
           'max_force': max_force,
           'rms_force': rms,
           'grand_total': gt,
           'grand_total_energy_per_atom': gtpa
           }
    return ret


def read_geometry_bas(filename):
    rf = open(filename, 'r')

    natom = int(rf.readline())
    symbols = []
    positions = np.zeros((natom, 3))

    for i in range(natom):
        line = rf.readline().split()
        symbols.append(atomic_symbol(int(line[0])))
        positions[i] = np.array(line[1:])

    return Structure(symbols=symbols, positions=positions, periodicity=False)


def write_geometry_bas(structure, filename):
    wf = open(filename, 'w')
    wf.write('%3d\n' % structure.natom)
    for i in range(structure.natom):
        if structure.is_periodic:
            x = structure.reduced[i, 0]
            y = structure.reduced[i, 1]
            z = structure.reduced[i, 2]
        else:
            x = structure.positions[i, 0]
            y = structure.positions[i, 1]
            z = structure.positions[i, 2]
        wf.write('%3d %13.7f %13.7f %13.7f\n' % (atomic_number(structure.symbols[i]), x, y, z))
    wf.close()


def get_fdata_info(fdata_path='Fdata'):
    rf = open(fdata_path + os.sep + 'info.dat')
    data = rf.read()
    ret = re.findall('([-. \d\w]*) - ([ \d\w;]*) \n', data)
    res = {}
    for i in ret:
        key = i[1]
        value = i[0]
        if key == 'Information for this species':
            cur_atom = int(value)
            res[cur_atom] = {}
        elif key == 'Element':
            res[cur_atom]['Element'] = value.strip()
        elif key == 'Nuclear Z':
            res[cur_atom]['Nuclear Z'] = int(value)
        elif key == 'Atomic Mass':
            res[cur_atom]['Atomic Mass'] = float(value)
        elif key == 'Number of shells; L for each shell':
            res[cur_atom]['Number of shells; L for each shell'] = int(value)
        elif key == 'Radial cutoffs PP':
            res[cur_atom]['Radial cutoffs PP'] = float(value)
        elif key == 'Atomic energy':
            res[cur_atom]['Atomic energy'] = float(value)
    return res


def read_fireball_in(fpath='fireball.in'):

    rf = open(fpath)
    data = rf.readlines()
    ret = {}

    for iline in data:
        if iline.startswith('&'):
            if iline.startswith('&END'):
                curkey = None
            else:
                curkey = iline[1:].strip()
        elif '=' in iline:
            varname = iline.split('=')[0].strip()
            value = iline.split('=')[1].strip()
            if curkey is not None and curkey not in ret:
                ret[curkey] = {}
            ret[curkey][varname] = convert_value(value)
    return ret


def read_eigen(fpath='eigen.dat'):

    rf = open(fpath)
    data = rf.read()
    # Number of Eigenvalues
    nval = int(data.split()[1])
    # The eigenvalues
    eigen = [float(x) for x in data.split()[7:]]
    assert(len(eigen) == nval)
    return eigen


def read_final_fireball_relax(fpath):

    output = read_fireball_stdout(fpath)
    ret = {'energetics': output['energetics'][-1],
           'forces': output['forces'][-1],
           'rms_force': output['rms_force'][-1],
           'max_force': output['max_force'][-1],
           'grand_total': output['grand_total'][-1],
           'grand_total_energy_per_atom': output['grand_total_energy_per_atom'][-1]}

    return ret


def convert_value(value):

    ret = None
    try:
        ret = int(value)
    except ValueError:
        try:
            ret = float(value)
        except ValueError:
            try:
                ret = [int(x) for x in value.split()]
            except ValueError:
                ret = value
    return ret


def read_param(fpath='param.dat'):
    rf = open(fpath)
    data = rf.readlines()
    ret = {}
    for line in data:
        if ':' in line:
            key = line.split(':')[0].strip()
            value = line.split(':')[1].strip()
            if value == '':
                curkey = key
            else:
                if curkey not in ret:
                    ret[curkey] = {}
                ret[curkey][key] = convert_value(value)
        else:
            if '=' not in line and len(line.strip()) > 0:
                curkey = line.strip()
    return ret


def read_lvs(fpath='input.lvs'):
    cell = np.zeros((3, 3))
    rf = open(fpath)
    data = rf.readlines()
    for i in range(len(data)):
        if len(data[i].split()) > 2:
            cell[i] = [float(x) for x in data[i].split()[:3]]
    return list(cell.flatten())
