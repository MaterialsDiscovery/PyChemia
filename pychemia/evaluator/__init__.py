import os as _os
import shutil as _shutil
import subprocess as _subprocess
from abc import ABCMeta, abstractmethod
from .function import FunctionEvaluator, FunctionObjectiveFunction
from .cluster import cluster_worker, cluster_evaluator, cluster_launcher
from .cluster_fireball import cluster_fb_evaluator, cluster_fb_launcher, cluster_fb_worker
from .direct_evaluator import DirectEvaluator
from .Fireball2PyChemiaDB import FireballCollector


def execute(basedir, command, script):
    """
    Utility that copy a given script and execute the given
    command inside the directory

    :param basedir: (str) Basedir of the directory with script
    :param command: (str) Command to execute
    :param script: (str) Script to call inside basedir
    """
    cwd = _os.getcwd()
    if not _os.path.isfile(basedir + '/' + _os.path.basename(script)):
        _shutil.copy(script, basedir)
    _os.chdir(basedir)
    print('Executing... ' + command + ' ' + script)
    _subprocess.call([command, script])
    _os.chdir(cwd)


class Evaluator:
    __metaclass__ = ABCMeta

    def __init__(self):
        self.population = None

    @abstractmethod
    def initialize(self, population):
        self.population = population

    @abstractmethod
    def run(self):
        pass

    @abstractmethod
    def stop(self):
        pass

    @property
    @abstractmethod
    def is_running(self):
        pass
