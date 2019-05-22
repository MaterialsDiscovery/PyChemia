[![Build Status](https://travis-ci.org/MaterialsDiscovery/PyChemia.svg?branch=master)](https://travis-ci.org/MaterialsDiscovery/PyChemia)
[![PyPI version](https://badge.fury.io/py/pychemia.svg)](https://badge.fury.io/py/pychemia)
[![Coverage Status](https://coveralls.io/repos/github/MaterialsDiscovery/PyChemia/badge.svg?branch=master)](https://coveralls.io/github/MaterialsDiscovery/PyChemia?branch=master)
[![Documentation Status](https://readthedocs.org/projects/pychemia/badge/?version=latest)](http://pychemia.readthedocs.io/en/latest/?badge=latest)

PyChemia
========

Python Materials Discovery Framework

![PyChemia](https://raw.githubusercontent.com/MaterialsDiscovery/PyChemia/master/docs/images/PyChemia_Small.png)

PyChemia is an open-source Python Library for materials structural search. The purpose of the initiative is to create
a method agnostic framework for materials discovery and design using a variety of methods from Minima Hoping to
Soft-computing based methods. PyChemia is also a library for data-mining, using several methods to discover interesting
candidates among the materials already processed.

The core of the library is the Structure python class, it is able to describe periodic and non-periodic structures.
As the focus of this library is structural search the class defines extensive capabilities to modify atomic structures.

The library includes capability to read and write in several ab-initio codes. At the level of DFT, PyChemia support
VASP, ABINIT and Octopus. At Tight-binding level development is in process to support DFTB+ and Fireball.
This allows the library to compute electronic-structure properties using state-of-the-art ab-initio software packages
and extract properties from those calculations.

PyChemia Requirements
=====================

Before installing PyChemia, you may need to first install a few critical dependencies

Mandatory
---------

1. Python >= 3.6
   The library is tested on Travis for versions of Python 3.6
   Support for Python 2.7 has been removed 

   https://travis-ci.org/MaterialsDiscovery/PyChemia

2. [Numpy](http://www.numpy.org/ "Numpy") >= 1.12
   NumPy is a fundamental package for any Python scientific library.
   Numpy arrays are essential for efficient array manipulation.
   Many distributions come with relatively old version of numpy,
   so you should install a more recent version.

3. [SciPy](http://scipy.org/ "SciPy") >= 0.18
   SciPy is used for many linear algebra and FFT calls
   Most distros comes with scipy 0.13 or below.
   PyChemia uses scipy.spatial that have been actively developed
   since version 0.15, we have tested PyChemia on 0.17.

4. [Spglib](http://spglib.sourceforge.net/) >= 1.9.9
   Spglib is used to determine symmetry groups for structures

5. [PyMongo](http://api.mongodb.org/python/current/) >= 3.4
   At least for structural search PyChemia relies strongly in MongoDB and its
   python driver. For the MongoDB server, any version beyond 3.0 should be fine.
   We have tested pychemia on MongoDB 3.2

   pip install pymongo

   or

   pip install pymongo --user


Optional Recommended
--------------------


1. [nose](https://nose.readthedocs.io/en/latest/) >= 1.3.7
   A python library for testing, simply go to the source directory and execute

   nosetests -v

2. [Matplotlib](http://matplotlib.org/  "Matplotlib") >= 1.2
   Used to plot band structures, densities of states and other 2D plots


Optional
--------

1. [Pandas](http://pandas.pydata.org/ "Pandas")
   Library for Data Analysis used by the datamining modules

2. [PyMC](http://pymc-devs.github.io/pymc/index.html)
   PyMC is a python module that implements Bayesian statistical models and fitting algorithms
   Important for the datamining capabilities of PyChemia

3. [Mayavi](http://docs.enthought.com/mayavi/mayavi/ "Mayavi") >= 4.1
   Some basic visualization tools are incorporated using this library

4. [ScientificPython](http://dirac.cnrs-orleans.fr/plone/software/scientificpython/overview/ "Scientific Python") >2.6
   This library is used for reading and writing NetCDF files

5. [pymatgen](http://www.pymatgen.org "pymatgen") >= 2.9
   pymatgen is an excellent library for materials analysis

6. [ASE](https://wiki.fysik.dtu.dk/ase/ "Atomic Simulation Environment")
   Atomic Simulation Environment is another good library for ab-initio calculations.
   Quite impressive for the number of ab-initio packages supported

7. [qmpy](http://oqmd.org/static/docs/index.html "qmpy")
   The Python library behind the [Open Quantum Materials Database](http://oqmd.org).
   The OQMD is a database of DFT calculated structures.
   For the time being the database contains more than 300000 structures, with more than
   90% of them with the electronic ground-state computed.

8. [coverage](https://bitbucket.org/ned/coveragepy) >= 4.0.1
    Provides code coverage analysis

9. [python-coveralls](https://github.com/z4r/python-coveralls)
    To submit coverage information to coveralls.io

    https://coveralls.io/github/MaterialsDiscovery/PyChemia

Documentation
=============

Instructions for installation, using and programming scripts with PyChemia
can be found on two repositories for documentation:

* Read The Docs:
   
   http://pychemia.readthedocs.io/en/latest
      
* Python Hosted:
    
   http://pythonhosted.org/pychemia

Structure of the Library
========================

![PyChemia](https://raw.githubusercontent.com/MaterialsDiscovery/PyChemia/master/docs/images/PyChemia_code.png)

![PyChemia](https://raw.githubusercontent.com/MaterialsDiscovery/PyChemia/master/docs/images/PyChemia_workflow.png)

Contributors
============

1. Prof. Aldo H. Romero [West Virginia University] (Project Director)

2. Guillermo Avendaño-Franco [West Virginia University] (Basic Infrastructure)

3. Adam Payne [West Virginia University] (Bug fixes (Populations, Relaxators, and KPoints) )

4. Irais Valencia Jaime [West Virginia University] (Simulation and testing)

5. Sobhit Singh [West Virginia University] (Data-mining)

6. Francisco Muñoz [Universidad de Chile] (PyPROCAR)

7. Wilfredo Ibarra Hernandez [West Virginia University] (Interface with MAISE)
