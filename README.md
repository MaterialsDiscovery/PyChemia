[![PyPI version](https://badge.fury.io/py/pychemia.svg)](https://badge.fury.io/py/pychemia)
[![Build Status](https://travis-ci.org/MaterialsDiscovery/PyChemia.svg?branch=master)](https://travis-ci.org/MaterialsDiscovery/PyChemia)
[![Documentation Status](https://readthedocs.org/projects/pychemia/badge/?version=latest)](http://pychemia.readthedocs.io/en/latest/?badge=latest)
[![Coverage Status](https://coveralls.io/repos/github/MaterialsDiscovery/PyChemia/badge.svg?branch=master)](https://coveralls.io/github/MaterialsDiscovery/PyChemia?branch=master)
[![HitCount](http://hits.dwyl.io/MaterialsDiscovery/PyChemia.svg)](http://hits.dwyl.io/MaterialsDiscovery/PyChemia)


PyChemia
========

Python Materials Discovery Framework

![PyChemia](https://raw.githubusercontent.com/MaterialsDiscovery/PyChemia/master/docs/images/PyChemia_Small.png)

PyChemia is an open-source Python Library for computational materials science. The library offers classes to manipulate crystal structures and compositions as well interfaces with first principles codes. Those abstractions allow the development of more complex tasks being built on top PyChemia basic building blocks.

PyChemia includes support to read and write inputs for several ab-initio codes. At the level of DFT, PyChemia supports VASP, ABINIT and Octopus. At Tight-binding level, PyChemia support DFTB+ and Fireball.
This allows the library to facilitate the calculations of electronic-structure properties using state-of-the-art ab-initio software packages.

When dealing with populations of structures or input parameters, PyChemia offers classes to stored and manipulate those collections. To achieve that PyChemia works by taking advantages of MongoDB, a document-based NoSQL database engine. 

PyChemia also implements several global optimization algorithms than when used together with an ab-initio code can perform a structural search of materials or some other applications where a local minimization is not enough.

Finally, Pychemia offers visualization and data mining capabilities useful to make sense from large amounts of calculations.

Installation
------------

Using pip:

    pip install pychemia --user

Using virtualenv

    virtualenv pychemia
    source pychemia/bin/activate
    pip install pychemia

You can also download the complete source code from GitHub

    git clone https://github.com/MaterialsDiscovery/PyChemia.git


Documentation
-------------

Documentation is hosted on [Read the Docs](https://readthedocs.org/projects/pychemia/) also available with Short URLs [readthedocs](http://pychemia.readthedocs.io) and [rtfd](http://pychemia.rtfd.io)

Documentation is also hosted on [Python Hosted](http://pythonhosted.org/pychemia/index.html)
