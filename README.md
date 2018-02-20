[![PyPI version](https://badge.fury.io/py/pychemia.svg)](https://badge.fury.io/py/pychemia)
[![Build Status](https://travis-ci.org/MaterialsDiscovery/PyChemia.svg?branch=master)](https://travis-ci.org/MaterialsDiscovery/PyChemia)
[![Documentation Status](https://readthedocs.org/projects/pychemia/badge/?version=latest)](http://pychemia.readthedocs.io/en/latest/?badge=latest)
[![Coverage Status](https://coveralls.io/repos/github/MaterialsDiscovery/PyChemia/badge.svg?branch=master)](https://coveralls.io/github/MaterialsDiscovery/PyChemia?branch=master)
[![HitCount](http://hits.dwyl.io/MaterialsDiscovery/PyChemia.svg)](http://hits.dwyl.io/MaterialsDiscovery/PyChemia)

PyChemia
========

Python Framework for Materials Discovery and Design

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

    $ pip install pychemia --user

Using virtualenv

    $ virtualenv venv_pcm
    $ source venv_pcm/bin/activate
    (venv_pcm)$ pip install pychemia

To deactivate the virtual environment execute:

    (venv_pcm)$ deactivate

For developers, or if you want to contribute with modifications, or bug corrections you can also download the complete source code from GitHub

    $ git clone https://github.com/MaterialsDiscovery/PyChemia.git

In that case you need to set your $PYTHONPATH before trying to import the library. For Korn Shell

    $ PYTHONPATH=<PATH_TO>/PyChemia:$PYTHONPATH
    $ export PYTHONPATH
    
On Bash (sh and bash) 

    $ export PYTHONPATH=<PATH_TO>/PyChemia:$PYTHONPATH

On c shell (csh and tcsh)

    $ setenv PYTHONPATH <PATH_TO>/PyChemia:$PYTHONPATH

Documentation
-------------

Documentation is hosted on [Read the Docs](https://readthedocs.org/projects/pychemia/) also available with Short URLs [readthedocs](http://pychemia.readthedocs.io) and [rtfd](http://pychemia.rtfd.io)

Documentation is also hosted on [Python Hosted](http://pythonhosted.org/pychemia/index.html)

Sources
-------

The main repository is on [GitHub](https://github.com/MaterialsDiscovery/PyChemia)

Sources and wheel binaries are also distrubuted on [PyPI](https://pypi.python.org/pypi/pychemia) or [PyPI](https://pypi.org/project/pychemia/)
