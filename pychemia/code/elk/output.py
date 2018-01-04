from ..codes import CodeOutput


class ElkOutput(CodeOutput):

    @property
    def is_finished(self):
        return False

    def read(self):
        pass

    def __init__(self):
        CodeOutput.__init__(self)
