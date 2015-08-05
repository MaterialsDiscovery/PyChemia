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
import shutil
import json
import tempfile
import subprocess

import pychemia
if pychemia.HAS_SCIPY and pychemia.HAS_SCIENTIFIC:
    import pychemia.code.abinit


def test_example2():
    """
    Example of a multiple calc          :
    """
    if not (pychemia.HAS_SCIPY and pychemia.HAS_SCIENTIFIC):
        return

    path = 'pychemia/test/data'
    assert(os.path.isdir(path))
    workdir = tempfile.mkdtemp()
    print(workdir)
    assert(os.path.isfile(path + '/abinit_04/t44.in'))
    var = pychemia.code.abinit.InputVariables(path + '/abinit_04/t44.in')
    abifiles = pychemia.code.abinit.AbiFiles(workdir)
    abifiles.set_input(var)
    abifiles.set_psps('LDA', 'FHI')
    abifiles.create()

    res = []
    wf = open(workdir + '/results.json', 'w')
    cwd = os.getcwd()
    os.chdir(workdir)

    for i in range(3):
        print(i)
        var.set_value('ecut', var.get_value('ecut') + 3)
        if i > 0:
            var.set_value('irdwfk', 1)
        var.write(abifiles.get_input_filename())
        abifile = open(workdir + '/abinit.files')
        logfile = open(workdir + '/abinit.log', 'w')
        if which('abinit') is not None:
            subprocess.call(['abinit'], stdin=abifile, stdout=logfile)
        else:
            print('The executable "abinit" is not in the PATH')
            print('Using the results of a previous calc')
        if os.path.isfile('abinit-o_WFK'):
            shutil.copyfile('abinit-o_WFK', 'abinit-i_WFK')
        data = pychemia.code.abinit.netcdf2dict(workdir + '/abinit-o_OUT.nc')
        os.rename(workdir + '/abinit-o_OUT.nc', '%s/abinit-o_OUT.nc_%d' % (workdir,i))
        res.append({'ecut': data['ecut'], 'etotal': data['etotal']})

    os.chdir(cwd)
    json.dump(res, wf)
    wf.close()

    if which('abinit') is None:
        res = json.load(open(path+os.sep+'abinit_04'+os.sep+'results.json'))

    assert (res[0]['etotal']+4.19954643154 < 1E-6)
    assert (res[1]['etotal']+4.19954643154 < 1E-6)
    assert (res[2]['etotal']+4.19954643154 < 1E-6)
    shutil.rmtree(workdir)


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
