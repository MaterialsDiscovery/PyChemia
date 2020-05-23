import os
import json
import subprocess
from setuptools import setup, find_packages, Extension
from distutils.command.sdist import sdist as _sdist
import pathlib

try:
    from Cython.Build import cythonize
    from Cython.Distutils import build_ext
except ImportError:
    USE_CYTHON = False
else:
    USE_CYTHON = True


# Return the git revision as a string
# Copied from scipy's setup.py
def git_version():
    def _minimal_ext_cmd(cmd):
        # construct minimal environment
        env = {}
        for k in ['SYSTEMROOT', 'PATH']:
            v = os.environ.get(k)
            if v is not None:
                env[k] = v
        # LANGUAGE is used on win32
        env['LANGUAGE'] = 'C'
        env['LANG'] = 'C'
        env['LC_ALL'] = 'C'
        out = subprocess.Popen(cmd, stdout=subprocess.PIPE, env=env).communicate()[0]
        return out

    try:
        out = _minimal_ext_cmd(['git', 'rev-parse', 'HEAD'])
        GIT_REVISION = out.strip().decode('ascii')
    except OSError:
        GIT_REVISION = "Unknown"

    return GIT_REVISION


def get_version_info():
    # Adding the git rev number needs to be done inside
    # write_version_py(), otherwise the import of scipy.version messes
    # up the build under Python 3.

    basepath=pathlib.Path(__file__).parent.absolute()
    print(basepath)

    rf = open(str(basepath)+os.sep+'setup.json')
    release_data = json.load(rf)
    rf.close()

    FULLVERSION = release_data['version']

    if os.path.exists('.git'):
        GIT_REVISION = git_version()
    elif os.path.exists('scipy/version.py'):
        # must be a source distribution, use existing version file
        # load it as a separate module to not load scipy/__init__.py
        import runpy
        ns = runpy.run_path('pychemia/version.py')
        GIT_REVISION = ns['git_revision']
    else:
        GIT_REVISION = "Unknown"

    if not ISRELEASED:
        FULLVERSION += '.dev0+' + GIT_REVISION[:7]

    return release_data, FULLVERSION, GIT_REVISION


def write_version_py(filename='pychemia/version.py'):
    cnt = """
# THIS FILE IS GENERATED FROM PYCHEMIA SETUP.PY
short_version = '%(version)s'
version = '%(version)s'
full_version = '%(full_version)s'
git_revision = '%(git_revision)s'
name = '%(name)s'
description = '%(description)s'
url = '%(url)s'
author = '%(author)s'
email = '%(email)s'
status = '%(status)s'
copyright = '%(copyright)s'
date = '%(date)s'
release = %(isrelease)s
if not release:
    version = full_version
"""
    release_data, FULLVERSION, GIT_REVISION = get_version_info()

    a = open(filename, 'w')
    try:
        a.write(cnt % {'version': release_data['version'],
                       'full_version': FULLVERSION,
                       'git_revision': GIT_REVISION,
                       'name': release_data['name'],
                       'description': release_data['description'],
                       'url': release_data['url'],
                       'author': release_data['author'],
                       'email': release_data['email'],
                       'status': release_data['status'],
                       'copyright': release_data['copyright'],
                       'date': release_data['date'],
                       'isrelease': str(ISRELEASED)})
    finally:
        a.close()
    return release_data


def get_scripts():
    return ['scripts' + os.sep + x for x in os.listdir('scripts') if x[-3:] == '.py']

###################################################################

ISRELEASED = False

KEYWORDS = ["electronic", "structure", "analysis", "materials", "discovery", "metaheuristics"]
CLASSIFIERS = [
    "Development Status :: 4 - Beta",
    "Intended Audience :: Science/Research",
    "License :: OSI Approved :: MIT License",
    "Natural Language :: English",
    "Operating System :: POSIX",
    "Programming Language :: Python",
    "Programming Language :: Python :: 3",
    "Programming Language :: Python :: 3.6",
    "Programming Language :: Python :: 3.7",
    "Topic :: Scientific/Engineering :: Chemistry",
    "Topic :: Scientific/Engineering :: Physics",
]
INSTALL_REQUIRES = ['numpy >= 1.17',
                    'scipy >= 1.3',
                    'spglib >= 1.9',
                    'pymongo >= 3.9',
                    'matplotlib >= 3.0',
                    'psutil >= 5.6']

###################################################################

print('Using Cython: %s' % USE_CYTHON)
data = write_version_py()

cmdclass = {}
ext = '.pyx' if USE_CYTHON else '.c'

class sdist(_sdist):
    def run(self):
            # Make sure the compiled Cython files in the distribution are
            # up-to-date
            from Cython.Build import cythonize
            cythonize(ext_modules, annotate=True, compiler_directives={'embedsignature': True})
            _sdist.run(self)

cmdclass['sdist'] = sdist

ext_modules = [Extension("pychemia.code.lennardjones.lj_utils", ['pychemia/code/lennardjones/lj_utils' + ext])]

if USE_CYTHON:
    cmdclass.update({'build_ext': build_ext})
    from Cython.Build import cythonize
    ext_modules = cythonize([Extension("pychemia.code.lennardjones.lj_utils", ['pychemia/code/lennardjones/lj_utils' + ext])])


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
    cmdclass=cmdclass,
    ext_modules=ext_modules
#    ext_modules=cythonize(ext_modules, annotate=True, compiler_directives={'embedsignature': True})
)

# ext_modules=cythonize(Extension('pychemia.code.lennardjones.hello',
#                                     ['pychemia/code/lennardjones/lj_utils.pyx'],
#                                     language='c',
#                                     extra_compile_args='-march=native'))
