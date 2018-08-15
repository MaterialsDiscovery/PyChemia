from ..codes import CodeRun
from .input import AbinitInput
import os


class AbinitRun(CodeRun):

    def __init__(self, workdir):
        CodeRun.__init__(self, executable='abinit', workdir=workdir, use_mpi=True)
        self.stdout_filename = 'abinit.out'
        self.stderr_filename = 'abinit.err'
        self.stdin_filename = 'abinit.files'

    def set_inputs(self):
        if not os.path.isfile(self.workdir+os.sep+'abinit.in'):
            raise ValueError("Input file 'abinit.in' not found at: %s" % self.workdir)
        else:
            self.input_path = self.workdir + os.sep + 'abinit.in'
            self.input = AbinitInput(self.input_path)

    def get_outputs(self):
        pass

    def is_finished(self):
        pass
    
    def is_successful(self):
        pass
