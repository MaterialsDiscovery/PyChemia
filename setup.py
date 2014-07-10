from distutils.core import setup

setup(
    name='pychemia',
    version='0.1.0',
    author='Guillermo Avendano Franco',
    author_email='gtux.gaf@gmail.com',
    packages=['pychemia',
              'pychemia.external.ase',
              'pychemia.external.pymatgen'],
    url='http://pypi.python.org/pypi/pychemia/',
    license='LICENSE.txt',
    description='Python framework for Materials Discovery and Design',
    long_description=open('README.md').read(),
    install_requires=["numpy >= 1.5",
                      "scipy >= 0.9",
                      "pymatgen >= 2.9",
                      "matplotlib >= 1.2",
                      "mayavi >= 4.1",
                      "ScientificPython >2.6",
                      "ase"],
)
