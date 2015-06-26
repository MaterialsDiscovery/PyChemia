__author__ = 'Guillermo Avendano-Franco'

import os
import subprocess
from pychemia import pcm_log
from pychemia.code import Codes
from pychemia.utils.periodic import atomic_number


class FireBall(Codes):
    def __init__(self):
        self.workdir = None
        # The five sections of Fireball
        self.option = {}
        self.tds = {}
        self.quench = {}
        self.tds = {}
        self.mesh = {}

        self.fdata_path = None
        self.structure = None
        self.binary = 'fireball.x'
        self.runner = None
        self.kpoints = None
        Codes.__init__(self)

    def initialize(self, workdir, structure, kpoints, binary='fireball.x'):
        assert structure.is_crystal
        assert structure.is_perfect
        self.structure = structure
        self.workdir = workdir
        if not os.path.lexists(workdir):
            os.mkdir(workdir)
        self.kpoints = kpoints
        self.binary = binary

    def set_inputs(self):
        self.write_input(filename=self.workdir + os.sep + 'fireball.in')
        self.write_basis(filename=self.workdir + os.sep + 'input.bas')
        self.write_lattice(filename=self.workdir + os.sep + 'input.lvs')
        self.write_kpoints(filename=self.workdir + os.sep + 'input.kpts')
        self.link_fdata()

    def get_outputs(self):
        pass

    def run(self, use_mpi=False, omp_max_threads=0, mpi_num_procs=1):
        cwd = os.getcwd()
        os.chdir(self.workdir)
        stdout = open('fireball.log', 'w')
        sp = subprocess.Popen(self.binary, stdout=stdout)
        os.chdir(cwd)
        self.runner = sp
        return sp

    def run_status(self):
        if self.runner is None:
            pcm_log.info('Fireball not finish')
            filename = self.workdir + os.sep + 'fireball.log'
            if os.path.exists(filename):
                results = read_fireball_stdout(filename=filename)
            return
        if self.runner.poll() == 0:
            pcm_log.info('Fireball complete normally')

    def finalize(self):
        pass

    def write_input(self, filename='fireball.in'):

        wf = open(filename, 'w')

        if self.option:
            wf.write('@OPTION\n')
            for variable in self.option:
                wf.write('%s = %s\n' % (variable, str(self.option[variable])))
            wf.write('@END\n')

        if self.output:
            wf.write('@OUTPUT\n')
            for variable in self.output:
                wf.write('%s = %s\n' % (variable, str(self.output[variable])))
            wf.write('@END\n')

        if self.quench:
            wf.write('@QUENCH\n')
            for variable in self.quench:
                wf.write('%s = %s\n' % (variable, str(self.quench[variable])))
            wf.write('@END\n')

        if self.mesh:
            wf.write('@MESH\n')
            for variable in self.mesh:
                wf.write('%s = %s\n' % (variable, str(self.mesh[variable])))
            wf.write('@END\n')

        if self.tds:
            wf.write('@TDS\n')
            for variable in self.tds:
                wf.write('%s = %s\n' % (variable, str(self.tds[variable])))
            wf.write('@END\n')

        wf.close()

    def write_basis(self, filename='input.bas'):

        wf = open(filename, 'w')
        wf.write('%d\n' % self.structure.natom)
        for i in range(self.structure.natom):
            x = self.structure.reduced[i, 0]
            y = self.structure.reduced[i, 1]
            z = self.structure.reduced[i, 2]
            wf.write('%d %13.7f %13.7f %13.7f\n' % (atomic_number(self.structure.symbols[i]), x, y, z))
        wf.close()

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
        linkpath = self.workdir + os.sep + 'Fdata'
        if os.path.exists(linkpath):
            os.remove(linkpath)
        os.symlink(self.fdata_path, linkpath)


def read_fireball_stdout(filename):
    rf = open(filename, 'r')
    return rf
