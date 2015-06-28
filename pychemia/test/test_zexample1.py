#!/usr/bin/env python
"""
This example shows how to automatize the execution
of ABINIT starting with just the input file "t44.in"

The abinit.files is written with the right
locations for the pseudopotentials and the results
are post-process using the output in NetCDF format
"abinit-o_OUT.nc"
"""

import os
import sys
import subprocess

import pychemia
if pychemia.HAS_SCIPY and pychemia.HAS_SCIENTIFIC:
    import pychemia.code.abinit as pa


path = 'pychemia/test/data'


def test_example1():
    """
    Example of a simple calc            :
    """
    if pychemia.HAS_SCIPY and pychemia.HAS_SCIENTIFIC:
        assert(os.path.isdir(path))
        workdir = get_path()
        print(workdir)
        creation(workdir)
        execution(workdir)
        ecut, etotal = datamining(workdir)
        print('ecut=', ecut)
        print('etotal=', etotal)
        assert ecut == 15.0
        assert abs(etotal + 4.19934820332) < 0.001


def get_path():
    workdir = None
    if os.path.isdir(path + '/abinit_04'):
        workdir = path + '/abinit_04'
    elif os.path.isdir('data/abinit_04'):
        workdir = 'data/abinit_04'
    else:
        print('The directory "abinit_04" was not found')
        exit(1)
    return workdir


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
        for ipath in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(ipath, program)
            if is_exe(exe_file):
                return exe_file
    return None


def creation(filep):
    """
    Create the input file from the original and
    set the proper values for abinit.files
    """
    var = pa.InputVariables(filep + '/t44.in')
    abifile = pa.AbiFiles(filep)
    abifile.set_input(var)
    abifile.set_psps('LDA', 'FHI')
    var.write(abifile.get_input_filename())
    print(abifile)
    print(var)
    abifile.create()


def execution(filep):
    """
    Execute ABINIT in the given path
    """
    abifile = open(filep + '/abinit.files')
    logfile = open(filep + '/abinit.log', 'w')
    cwd = os.getcwd()
    os.chdir(filep)
    if which('abinit') is not None:
        subprocess.call(['abinit'], stdin=abifile, stdout=logfile)
    else:
        print('The executable "abinit" is not in the PATH')
        print('Using the results of a previous calc')
    files = [x for x in os.listdir('.') if x not in ['t44.in', 'abinit-o_OUT.nc']]
    for i in files:
        os.remove(i)
    os.chdir(cwd)
    abifile.close()
    logfile.close()


def datamining(filep):
    """
    Read some output variables from abinit-o_OUT.nc
    """
    print('Extracting ecut and etotal')
    data = pa.netcdf2dict(filep + '/abinit-o_OUT.nc')
    print(data)
    print(data['ecut'][0], data['etotal'][0])
    return data['ecut'][0], data['etotal'][0]


if __name__ == '__main__':
    path = get_path()
    if len(sys.argv) == 1:
        test_example1()
    elif len(sys.argv) == 2:
        if sys.argv[1] == 'create':
            creation(path)
        if sys.argv[1] == 'execute':
            execution(path)
        elif sys.argv[1] == 'result':
            datamining(path)
