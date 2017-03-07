from setuptools import setup, find_packages

setup(
    name='pychemia',
    version='0.17.3',
    author='Guillermo Avendano-Franco',
    author_email='gufranco@mail.wvu.edu',
    packages=find_packages(),
    url='https://github.com/MaterialsDiscovery/PyChemia',
    license='LICENSE.txt',
    description='Python framework for Materials Discovery and Design',
    long_description=open('README.md').read(),
    install_requires=['numpy >= 1.12.0',
                      'scipy >= 0.18.0',
                      'future>=0.16.0',
                      'spglib>=1.9.9',
                      'pymongo>=3.4.0'],
    keywords=["electronic", "structure", "analysis", "materials", "discovery", "metaheuristics"],
)
