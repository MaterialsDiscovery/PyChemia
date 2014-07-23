import os as _os
import shutil as _shutil
import subprocess as _subprocess


def execute(basedir, command, script):
    """
    Utility that copy a given script and execute the given
    command inside the directory
    """
    cwd = _os.getcwd()
    if not _os.path.isfile(basedir + '/' + _os.path.basename(script)):
        _shutil.copy(script, basedir)
    _os.chdir(basedir)
    print('Executing... ' + command + ' ' + script)
    _subprocess.call([command, script])
    _os.chdir(cwd)
