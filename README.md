PyChemia
========

Python Materials Discovery Framework

<center>
![PyChemia](https://raw.githubusercontent.com/MaterialsDiscovery/PyChemia/master/docs/images/PyChemia_Small.png)
</center>

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

PyChemia requirements
---------------------

Before installing PyChemia, you may need to first install a few critical dependencies

### Mandatory

1. Python >= 2.7.x or Python >= 3.1.x
   The library is tested for the current version of Python

2. [Numpy](http://www.numpy.org/ "Numpy") >= 1.5
   NumPy is a fundamental package for any Python scientific library.
   Numpy arrays are essential for efficient array manipulation.

3. [SciPy](http://scipy.org/ "SciPy") >= 0.9
   SciPy is used for many linear algebra and FFT calls

4. [spglib](http://spglib.sourceforge.net/)
   Spglib is used to determine symmetry groups for structures

5. [pymongo](http://api.mongodb.org/python/current/)
   At least for structural search PyChemia relies strongly in MongoDB and its
   python driver

### Optional

1. [Pandas](http://pandas.pydata.org/ "Pandas")
   Library for Data Analysis used by the datamining modules

1. [PyMC](http://pymc-devs.github.io/pymc/index.html)
   PyMC is a python module that implements Bayesian statistical models and fitting algorithms
   Important for the datamining capabilities of PyChemia

1. [Matplotlib](http://matplotlib.org/  "Matplotlib") >= 1.2
   Used to plot band structures, densities of states and other 2D plots

1. [Mayavi](http://docs.enthought.com/mayavi/mayavi/ "Mayavi") >= 4.1
   Some basic visualization tools are incorporated using this library

1. [ScientificPython](http://dirac.cnrs-orleans.fr/plone/software/scientificpython/overview/ "Scientific Python") >2.6
   This library is used for reading and writing NetCDF files

1. [pymatgen](http://www.pymatgen.org "pymatgen") >= 2.9
   pymatgen is an excellent library for materials analysis

1. [ASE](https://wiki.fysik.dtu.dk/ase/ "Atomic Simulation Environment")
   Atomic Simulation Environment is another good library for ab-initio calculations.
   Quite impressive for the number of ab-initio packages supported

1. [qmpy](http://oqmd.org/static/docs/index.html "qmpy")
   The Python library behind the [Open Quantum Materials Database](http://oqmd.org).
   The OQMD is a database of DFT calculated structures.
   For the time being the database contains more than 300000 structures, with more than
   90% of them with the electronic ground-state computed.

How to install PyChemia
-----------------------

The best way to install PyChemia is from its GitHub repository

    git clone https://github.com/MaterialsDiscovery/PyChemia.git

Once the repository is downloaded, you can install it using

    cd PyChemia
    python setup.py built
    sudo python setup.py install

Structure of the Library
------------------------

![PyChemia](https://raw.githubusercontent.com/MaterialsDiscovery/PyChemia/master/docs/images/PyChemia_code.png)

![PyChemia](https://raw.githubusercontent.com/MaterialsDiscovery/PyChemia/master/docs/images/PyChemia_workflow.png)

Contributors
------------

1. Prof. Aldo H. Romero [West Virginia University] (Project Director)

1. Guillermo Avendaño-Franco [West Virginia University] (Basic Infraestructure)

1. Irais Valencia Jaime [West Virginia University] (Simulation and testing)

1. Sobhit Singh [West Virginia University] (Data-mining)

1. Francisco Muñoz (Universidad de Chile)