import os
from setuptools import setup, find_packages, Extension
import json

try:
    from Cython.Build import cythonize
    from Cython.Distutils import build_ext
except ImportError:
    USE_CYTHON = False
else:
    USE_CYTHON = True


rf = open('pychemia' + os.sep + 'setup.json')
data = json.load(rf)
rf.close()


def get_scripts():
    return ['scripts' + os.sep + x for x in os.listdir('scripts') if x[-3:] == '.py']


###################################################################


KEYWORDS = ["electronic", "structure", "analysis", "materials", "discovery", "metaheuristics"]
CLASSIFIERS = [
    "Development Status :: 4 - Beta",
    "Intended Audience :: Science/Research",
    "License :: OSI Approved :: MIT License",
    "Natural Language :: English",
    "Operating System :: POSIX",
    "Programming Language :: Python",
    "Programming Language :: Python :: 3",
    "Programming Language :: Python :: 3.5",
    "Programming Language :: Python :: 3.6",
    "Programming Language :: Python :: 3.7",
    "Topic :: Scientific/Engineering :: Chemistry",
    "Topic :: Scientific/Engineering :: Physics",
]
INSTALL_REQUIRES = ['numpy >= 1.12.0',
                    'scipy >= 0.18.0',
                    'spglib >= 1.9.9',
                    'pymongo >= 3.4.0',
                    'psutil >= 5.4.7']

###################################################################

ext = '.pyx' if USE_CYTHON else '.c'

extensions = [Extension("pychemia.code.lennardjones.lj_utils", ['pychemia/code/lennardjones/lj_utils' + ext])]

setup(
    name=data['name'],
    version=data['version'],
    author=data['author'],
    author_email=data['email'],
    packages=find_packages(exclude=['scripts', 'docs', 'tests']),
    url=data['url'],
    license='LICENSE.txt',
    description=data['description'],
    long_description=open('README').read(),
    install_requires=INSTALL_REQUIRES,
    keywords=KEYWORDS,
    classifiers=CLASSIFIERS,
    python_requires='>=2.7, !=3.0.*, !=3.1.*, !=3.2.*, <4',
    package_data={'': ['setup.json']},
    scripts=get_scripts(),
    ext_modules=cythonize(extensions, annotate=True, compiler_directives={'embedsignature': True})
)

# ext_modules=cythonize(Extension('pychemia.code.lennardjones.hello',
#                                     ['pychemia/code/lennardjones/lj_utils.pyx'],
#                                     language='c',
#                                     extra_compile_args='-march=native'))
