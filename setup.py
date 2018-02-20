import os
from setuptools import setup, find_packages


def get_scripts():
    return ['scripts' + os.sep + x for x in os.listdir('scripts') if x[-3:] == '.py']


###################################################################

NAME = "pychemia"
VERSION = '0.18.2.20'
KEYWORDS = ["electronic", "structure", "analysis", "materials", "discovery", "metaheuristics"]
CLASSIFIERS = [
    "Development Status :: 4 - Beta",
    "Intended Audience :: Science/Research",
    "License :: OSI Approved :: MIT License",
    "Natural Language :: English",
    "Operating System :: POSIX",
    "Programming Language :: Python",
    "Programming Language :: Python :: 2",
    "Programming Language :: Python :: 2.7",
    "Programming Language :: Python :: 3",
    "Programming Language :: Python :: 3.4",
    "Programming Language :: Python :: 3.5",
    "Programming Language :: Python :: 3.6",
    "Topic :: Scientific/Engineering :: Chemistry",
    "Topic :: Scientific/Engineering :: Physics",
]
INSTALL_REQUIRES = ['numpy >= 1.12.0',
                    'scipy >= 0.18.0',
                    'future>=0.16.0',
                    'spglib>=1.9.9',
                    'pymongo>=3.4.0']

###################################################################

setup(
    name=NAME,
    version=VERSION,
    author='Guillermo Avendano-Franco',
    author_email='gufranco@mail.wvu.edu',
    packages=find_packages(exclude=['scripts', 'docs', 'tests']),
    url='https://github.com/MaterialsDiscovery/PyChemia',
    license='LICENSE.txt',
    description='Python framework for Materials Discovery and Design',
    long_description=open('README').read(),
    install_requires=INSTALL_REQUIRES,
    keywords=KEYWORDS,
    classifiers=CLASSIFIERS,
    python_requires='>=2.7, !=3.0.*, !=3.1.*, !=3.2.*, <4',
    data_files=[('my_data', ['tests/data/Au.cif'])],
    scripts=get_scripts(),
)


