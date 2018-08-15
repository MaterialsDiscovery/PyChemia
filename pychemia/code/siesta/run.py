from ..codes import CodeRun
from .input import SiestaInput
import os


class SiestaRun(CodeRun):

    def __init__(self, workdir, input_path, pseudo_path, pseudo_file=None, pseudo_list=None):
        CodeRun.__init__(self, executable='siesta', workdir=workdir, use_mpi=False)

        if pseudo_file is None and pseudo_list is None:
            raise ValueError("Declare either the list of pseudo names or a file that contain the list.")

        if pseudo_list is not None and not os.path.isdir(pseudo_path):
            raise ValueError('ERROR: Directory with Pseudo-potentials not found: %s' % pseudo_path)

        if pseudo_path is not None and not os.path.exists(input_path):
            raise ValueError('ERROR: File not found: %s' % input_path)

        self.workdir = workdir
        self.input_path = input_path
        self.stdin_filename = 'siesta.fdf'
        self.stdout_filename = 'siesta.out'
        self.stderr_filename = 'siesta.err'

        self.pseudo_path = os.path.abspath(pseudo_path)

        if os.path.exists(pseudo_file):
            rf = open(pseudo_file)
            self.pseudos = rf.read().split()

    def set_inputs(self):
        # Create a siesta input object
        si = SiestaInput(self.input_path)
        
        if not os.path.exists(self.workdir):
            os.mkdir(self.workdir)

        si.write(self.workdir+os.sep+'siesta.fdf')

        for i in self.pseudos:
            psf = self.workdir + os.sep + i + '.psf'
            if os.path.isfile(psf):
                os.remove(psf)
            os.symlink(self.pseudo_path + os.sep + i + '.psf', psf)

    def get_outputs(self):
        pass

    def is_finished(self):
        outputfile = self.workdir+os.sep+'siesta.out'
        if not os.path.isfile(outputfile):
            return False
        rf = open(outputfile)
        data = rf.read()
        if data[-4:] == 'Job completed\n':
            return True
        else:
            return False

    def is_successful(self):
        return self.is_finished()

