import os
from .abifiles import AbiFiles
from .input import AbinitInput
from ..codes import CodeRun
from pychemia import pcm_log


class AbinitJob(CodeRun):
    def __init__(self, executable='abinit', workdir='.'):
        CodeRun.__init__(self, executable=executable, workdir=workdir, use_mpi=True)
        self.abifile = None
        self.kpoints = None
        self.structure = None
        self.kind = None
        self.exchange = None
        self.inp = AbinitInput()
        self.stdout_filename = 'abinit.log'
        self.stdin_filename = 'abinit.files'
        self.stderr_filename = 'abinit.err'
        pcm_log.debug("Created AbinitJob with executable=%s and workdir=%s" % (executable, workdir))

    def finalize(self):
        self.stdout_file.close()

    def get_outputs(self):
        pass

    def initialize(self, structure=None, input_file='abinit.in', psp_kind='PBE', psp_exchange='ONC'):
        """
        Initialize the mandatory variables for a AbinitJob.
        The structure, input and pseudopotentials can change for the same AbinitJob.

        :param psp_kind: (str) Source of Pseudopotentials
        :param psp_exchange: (str) 'LDA' or 'GGA'
        :param structure: (pychemia.Structure) A pychemia structure for the input
        :param input_file: (str) Input file for ABINIT

        :return:
        """
        self.abifile = AbiFiles(self.workdir)
        if not os.path.lexists(self.workdir):
            os.mkdir(self.workdir)
        if structure is not None:
            pcm_log.debug("Initializing AbinitJob with structure=%s" % (structure.formula))
            assert structure.is_crystal
            assert structure.is_perfect
            self.structure = structure
            self.inp.from_structure(self.structure)
        if os.path.isfile(input_file):
            self.inp = AbinitInput(input_file)
        self.abifile.set_input(self.inp)
        self.kind = psp_kind
        self.exchange = psp_exchange

    def job_ion_relax(self, internal_opt=True, external_opt=True, ionmov=2, nstep=20, ntime=30,
                      tolmxf=1E-7, tolrff=1E-3, dilatmx=1.05, chkprim=0):

        pcm_log.debug("Setting variables for relaxation tolmxf=%e tolrff=%e" % (tolmxf, tolrff))

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
        self.inp.set_value('chkprim', chkprim)

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

    def set_psps(self, exchange='PBE', kind='ONC'):
        self.abifile.set_psps(exchange=exchange, kind=kind)

    def set_inputs(self):
        """
        Prepare all the inputs before running ABINIT

        :return:
        """
        self.set_psps(kind=self.kind, exchange=self.exchange)
        self.set_kpoints()
        self.write_abifiles()
        self.write_abiinput()

    def set_kpoints(self, kpoints=None):
        self.kpoints = kpoints
        if 'ngkpt' not in self.inp.variables.keys():
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
        self.inp.write(self.workdir + os.sep + 'abinit.in')

    def write_abifiles(self):
        self.abifile.create()

    def write_all(self):
        self.write_abifiles()
        self.write_abiinput()

