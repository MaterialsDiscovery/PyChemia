__author__ = 'Guillermo Avendano Franco'

import os

from pychemia.code._codes import Codes


class FireBall(Codes):

    def __init__(self):
        self._dirpath = None

    def initialize(self, dirpath):
        self._dirpath = dirpath
        if not os.path.lexists(dirpath):
            os.mkdir(dirpath)

    def set_inputs(self):
        pass

    def get_ouputs(self):
        pass

    def run(self):
        pass

    def finalize(self):
        pass

    @property
    def dirpath(self):
        return self._dirpath
