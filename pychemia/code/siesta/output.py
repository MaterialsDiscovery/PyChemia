from ..codes import CodeOutput


class SiestaOutput(CodeOutput):

    @property
    def is_finished(self):
        return False

    def __init__(self):
        CodeOutput.__init__(self)

    def read(self):
        pass
