import os
import subprocess
import pychemia

path = 'pychemia/test/data'


def notest_abinit2xyz():
    """
    Test command abinit2xyz                                      :
    """
    pychemia.info()
    print('Running abinit2xyz')
    script = 'bin/abi_abinit2xyz.py'
    args = path + '/abinit_03/abinit.in'
    command_line = "%s %s %s" % ('python', script, args)
    print(command_line)
    subprocess.check_output(command_line, shell=True)
    mol1 = pychemia.io.xyz.load(path + '/abinit_03/abinit_DS1.xyz')
    mol2 = pychemia.io.xyz.load(path + '/abinit_03/abinit_DS2.xyz')
    assert (mol1.natom == 16)
    assert (mol2.natom == 16)
    os.remove(path + '/abinit_03/abinit_DS1.xyz')
    os.remove(path + '/abinit_03/abinit_DS2.xyz')
    os.remove(path + '/abinit_03/abinit.files')


def notest_xyz2abinit():
    """
    Test command xyz2abinit                                      :
    """
    script = 'bin/abi_xyz2abinit.py'
    args = path + '/xyz/chlorophyll.xyz'
    command_line = "%s %s %s" % ('python', script, args)
    print(command_line)
    subprocess.check_output(command_line, shell=True)
    inp = pychemia.code.abinit.AbinitInput(path + '/xyz/chlorophyll.xyz.in')
    assert (inp.get_value('natom') == 140)
    os.remove(path + '/xyz/chlorophyll.xyz.in')
    return inp


def notest_plot_bonds():
    """
    Test command plot_bonds                                      :
    """
    script = 'bin/abi_plot_bonds.py'
    arg1 = path + '/abinit_01/abinit.files:11'
    arg2 = path + '/abinit_01/abinit.files:21'
    command_line = "%s %s %s %s" % ('python', script, arg1, arg2)
    print(command_line)
    subprocess.check_output(command_line, shell=True)
    os.remove("bonds.txt")
    os.remove("bonds.pdf")


def notest_plot_hist():
    """
    Test command plot_hist                                       :
    """
    script = 'bin/abi_plot_hist.py'
    arg1 = path + '/abinit_01/abinit.files'
    arg2 = '11'
    command_line = "%s %s %s %s" % ('python', script, arg1, arg2)
    print(command_line)
    subprocess.check_output(command_line, shell=True)
    os.remove(path + '/abinit_01/abinit-o_DS11.pdf')
