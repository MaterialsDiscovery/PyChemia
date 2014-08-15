import os
import subprocess

import pychemia

path = 'pychemia/test/data'

def test_abinit2xyz():
    """
    Test command abinit2xyz             :
    """
    script = 'bin/abi_abinit2xyz.py'
    args = path + '/abinit_03/abinit.in'
    subprocess.call(['python', script, args])
    mol1 = pychemia.io.xyz.load(path + '/abinit_03/abinit_DS1.xyz')
    mol2 = pychemia.io.xyz.load(path + '/abinit_03/abinit_DS2.xyz')
    assert (mol1.natom == 16)
    assert (mol2.natom == 16)
    os.remove(path + '/abinit_03/abinit_DS1.xyz')
    os.remove(path + '/abinit_03/abinit_DS2.xyz')
    os.remove(path + '/abinit_03/abinit.files')


def test_xyz2abinit():
    """
    Test command xyz2abinit             :
    """
    script = 'bin/abi_xyz2abinit.py'
    args = path + '/xyz/chlorophyll.xyz'
    subprocess.call(['python', script, args])
    inp = pychemia.code.abinit.InputVariables(path + '/xyz/chlorophyll.xyz.in')
    assert (inp.get_value('natom') == 140)
    os.remove(path + '/xyz/chlorophyll.xyz.in')
    return inp


def test_plot_bonds():
    """
    Test command plot_bonds             :
    """
    script = 'bin/abi_plot_bonds.py'
    arg1 = path + '/abinit_01/abinit.files:11'
    arg2 = path + '/abinit_01/abinit.files:21'
    subprocess.call(['python', script, arg1, arg2])
    os.remove("bonds.txt")
    os.remove("bonds.pdf")


def test_plot_hist():
    """
    Test command plot_hist              :
    """
    script = 'bin/abi_plot_hist.py'
    arg1 = path + '/abinit_01/abinit.files'
    arg2 = '11'
    subprocess.call(['python', script, arg1, arg2])
    os.remove(path + '/abinit_01/abinit-o_DS11.pdf')
