from setuptools import setup, find_packages, Extension

setup(
    name='pychemia',
    version='0.1.1',
    author='Guillermo Avendano Franco',
    author_email='gufranco@mail.wvu.edu',
    packages=find_packages(),
    url='https://github.com/MaterialsDiscovery/PyChemia',
    license='LICENSE.txt',
    description='Python framework for Materials Discovery and Design',
    long_description=open('README.md').read(),
    install_requires=["numpy >= 1.5", 'scipy', 'matplotlib'],
    keywords=["electronic", "structure", "analysis", "materials", "discovery", "metaheuristics"],
)
