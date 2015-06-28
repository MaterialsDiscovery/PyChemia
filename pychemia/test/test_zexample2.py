"""
This example shows how to automatize the execution
of ABINIT starting with just the input file "t44.in"

The abinit.files is written with the right
locations for the pseudopotentials and the results
are post-process using the output in NetCDF format
"abinit-o_OUT.nc"
"""

import os
import shutil
import subprocess

import pychemia
if pychemia.HAS_SCIPY and pychemia.HAS_SCIENTIFIC:
    import pychemia.code.abinit as pa
    import pychemia.code.abinit


def test_example2():
    """
    Example of a multiple calc          :
    """
    if not (pychemia.HAS_SCIPY and pychemia.HAS_SCIENTIFIC):
        exit(1)

    path = 'pychemia/test/data'
    assert(os.path.isdir(path))
    path = get_path()
    print(path)
    assert(os.path.isfile(path + '/t44.in'))
    var = pa.InputVariables(path + '/t44.in')
    abifiles = pa.AbiFiles(path)
    abifiles.set_input(var)
    abifiles.set_psps('LDA', 'FHI')
    abifiles.create()

    wf = open(path + '/results.txt', 'w')

    for i in range(3):
        print(i)
        var.set_value('ecut', var.get_value('ecut') + 3)
        if i > 0:
            var.set_value('irdwfk', 1)
        var.write(abifiles.get_input_filename())
        abifile = open(path + '/abinit.files')
        logfile = open(path + '/abinit.log', 'w')
        cwd = os.getcwd()
        os.chdir(path)
        if which('abinit') is not None:
            subprocess.call(['abinit'], stdin=abifile, stdout=logfile)
        else:
            print('The executable "abinit" is not in the PATH')
            print('Using the results of a previous calc')
        if os.path.isfile('abinit-o_WFK'):
            shutil.copyfile('abinit-o_WFK', 'abinit-i_WFK')
        os.chdir(cwd)
        print(os.getcwd())
        data = pychemia.code.abinit.netcdf2dict(path + '/abinit-o_OUT.nc')
        wf.write(str(data['ecut'][0]).rjust(20) + ' ' + str(data['etotal'][0]).rjust(20) + '\n')
        if os.path.isfile(path + '/abinit.out'):
            os.remove(path + '/abinit.out')
        if os.path.isfile(path + '/abinit.log'):
            os.remove(path + '/abinit.log')

    wf.close()
    if os.path.isfile(path + '/abinit.files'):
        os.remove(path + '/abinit.files')


def get_path():
    path = None
    if os.path.isdir('pychemia/test/data/abinit_04'):
        path = 'pychemia/test/data/abinit_04'
    elif os.path.isdir('data/abinit_04'):
        path = 'data/abinit_04'
    else:
        print('The directory "abinit_04" was not found')
        exit(1)
    return path


def which(program):
    """
    Search for the presence of an executable
    Found in: http://stackoverflow.com/questions/377017/test-if-executable-exists-in-python
    """
    import os

    def is_exe(filep):
        return os.path.isfile(filep) and os.access(filep, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file
    return None


if __name__ == '__main__':
    test_example2()
