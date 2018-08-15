from ..codes import CodeRun


class ElkRun(CodeRun):

    def __init__(self, workdir):
        CodeRun.__init__(self, executable='elk', workdir=workdir, use_mpi=True)

    def set_inputs(self):
        pass

    def get_outputs(self):
        pass

    def is_finished(self):
        pass

    def is_successful(self):
        pass
