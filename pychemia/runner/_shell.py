import os
import shutil
import socket
import subprocess
import xml.etree.ElementTree as ElementTree
from pathlib import Path

class SHELL_Runner:
    """
    Class to create and launch submission scripts for Torque/PBS

    """
    def __init__(self, workdir=None, submission_script='batch.slurm', slurm_params=None, input_files=None, jobname=None, command_line=None, use_symlinks=False):
        """
        """
        pass

