"""
Visualization of Anthracene
"""

import os
import time
import pychemia
import pychemia.external.pymatgen

path = 'tests/data'


def test_anthracene():
    """
    Visualization of anthracene         :
    """
    print(os.getcwd())
    assert (os.path.isdir(path))
    filename = path + '/xyz/anthracene.xyz'
    assert (os.path.isfile(filename))
    mol = pychemia.io.xyz.load(filename)
    assert (mol.natom == 24)
    assert (mol.is_crystal is False)
    assert (mol.is_periodic is False)
    fig = mol.plot()
    time.sleep(5)
    fig.remove()
    fig.stop()


def test_chlorophyll():
    """
    Visualization of chlorophyll        :
    """
    print(os.getcwd())
    assert (os.path.isdir(path))
    filename = path + '/xyz/chlorophyll.xyz'
    assert (os.path.isfile(filename))
    mol = pychemia.io.xyz.load(filename)
    assert (mol.natom == 140)
    assert (mol.is_periodic is False)
    assert (mol.is_crystal is False)
    fig = mol.plot()
    time.sleep(5)
    fig.remove()
    fig.stop()


def test_gold():
    """
    Visualization of Gold FCC lattice   :
    """
    a = 4.05
    b = a / 2
    fcc = pychemia.Structure(symbols=['Au'], cell=[[0, b, b], [b, 0, b], [b, b, 0]], periodicity=True)
    assert (fcc.natom == 1)
    assert (fcc.is_periodic is True)
    assert (fcc.is_crystal is True)
    fig = fcc.plot()
    time.sleep(5)
    fig.remove()
    fig.stop()
