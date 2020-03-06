from ..codes import CodeRun
from .input import OctopusInput
import os


class OctopusRun(CodeRun):

    def __init__(self, workdir='.'):
        CodeRun.__init__(self, executable='octopus', workdir=workdir, use_mpi=True)
        self.stdout_filename = 'octopus.out'
        self.stderr_filename = 'octopus.err'

    def set_inputs(self):
        if not os.path.isfile(self.workdir+os.sep+'inp'):
            raise ValueError("Input file 'inp' not found at: %s" % self.workdir)
        else:
            self.input_path = self.workdir + os.sep + 'inp'
            self.input = OctopusInput(self.input_path)

    def get_outputs(self):
        pass

    def is_finished(self):
        pass
    
    def is_successful(self):
        pass
