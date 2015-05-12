__author__ = 'Guillermo Avendano Franco'

import os

from _abifiles import AbiFiles
from _input import InputVariables
from pychemia.code import Codes


class AbinitJob(Codes):

    def __init__(self):
        Codes.__init__(self)
        self.abifile = AbiFiles(self.workdir)
        self.kpoints = None
        self.structure = None
        self.workdir = None
        self.inp = InputVariables()
        self.stdout_filename = 'abinit.log'
        self.stdin_filename = 'abinit.files'
        self.stderr_filename = 'abinit.err'

    def finalize(self):
        self.stdout_file.close()

    def get_outputs(self):
        pass

    def initialize(self, workdir, structure, kpoints, binary='abinit'):
        """
        Initialize the mandatory variables for a AbinitJob

        :param workdir: (str) The directory where the input files will be
        :param structure: (pychemia.Structure) A pychemia structure for the input
        :param kpoints: (pychemia.dft.KPoints) The Kpoints object to use
        :param binary: (str) The name or path for the ABINIT binary file

        :return:
        """
        assert structure.is_crystal
        assert structure.is_perfect
        self.structure = structure
        self.workdir = workdir
        if not os.path.lexists(workdir):
            os.mkdir(workdir)
        self.binary = binary
        self.kpoints = None

        self.inp.from_structure(self.structure)
        self.abifile.set_input(self.inp)

    def job_ion_relaxation(self, internal_opt=True, external_opt=True, ionmov=2, nstep=20, ntime=30,
                           tolmxf=1E-7, tolrff=1E-3, dilatmx=1.05):

        self.inp.clean()
        self.inp.from_structure(self.structure)

        if internal_opt and external_opt:
            self.inp.set_value('ndtset', 2)
            self.inp.set_value('optcell', 0, idtset=1)
            self.inp.set_value('optcell', 2, idtset=2)
            self.inp.set_value('ecutsm', 0.5, idtset=2)
        elif external_opt:
            self.inp.set_value('optcell', 2)

        self.inp.set_value('ionmov', ionmov)
        self.inp.set_value('nstep', nstep)
        self.inp.set_value('ntime', ntime)
        self.inp.set_value('tolmxf', tolmxf)
        self.inp.set_value('tolrff', tolrff)
        self.inp.set_value('dilatmx', dilatmx)

    def job_static(self):
        """
        Prepares a static calculation.

        :return:
        """
        self.inp.clean()
        self.inp.from_structure(self.structure)
        self.inp.set_value('ionmov', 0)
        self.inp.set_value('toldfe', 1E-7)

    def job_ecut_convergence(self):

        self.inp.clean()
        self.job_static()
        self.inp.set_value('ndtset', 8)
        self.inp.set_value('ecut1', 50)
        self.inp.set_value('ecut2', 100)
        self.inp.set_value('ecut3', 150)
        self.inp.set_value('ecut4', 200)
        self.inp.set_value('ecut5', 300)
        self.inp.set_value('ecut6', 500)

    def set_psps(self, exchange='LDA', kind='FHI'):
        self.abifile.set_psps(exchange=exchange, kind=kind)

    def set_inputs(self):
        """
        Prepare all the inputs before running ABINIT

        :return:
        """
        self.set_psps()
        self.set_kpoints()
        self.write_abifiles()
        self.write_abiinput()

    def set_kpoints(self):
        if self.kpoints is None:
            self.inp.set_value('kptrlen', 1)
        elif self.kpoints.kmode != 'gamma':
            raise ValueError('Not Implemented yet')
        else:
            self.inp.set_value('ngkpt', list(self.kpoints.grid))
            if self.kpoints.shifts is not None:
                self.inp.set_value('nshiftk', len(self.kpoints.shifts))
                self.inp.set_value('shiftk', self.kpoints.shifts)

    def set_ecut(self, ecut=50):
        self.inp.set_value('ecut', ecut)

    def write_abiinput(self):
        self.inp.write(self.workdir+os.sep+'abinit.in')

    def write_abifiles(self):
        self.abifile.create()
