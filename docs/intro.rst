Introduction
============

PyChemia is a open-source Python library for High-throughput first-principles
materials discovery. The focus of this library is on structural search and 
data analysis.
The ultimate purpose of the code is to optimize the search of new materials 
using a variety of methods such as Minima hopping method (MHM), soft-computing-based methods and statistical methods.

The main objectives of the code are:

1. Provide flexible classes for atomic structures such as molecules, clusters,
   thin films and crystals.

2. Manipulate both input and output for DFT and Tight-binding codes such as
   VASP. ABINIT, Fireball and DFTB+

3. Offer a robust architecture for storing a large collection of structures.
   Structural search methods generate many structures that are stored in
   repositories. Calculations done on those structures are also store in
   repositories.

4. Similarity analysis based on fingerprints, pair correlation and comparators.

5. Stability analysis for crystals. Including thermal analysis (Enthalpy) and
   dynamic stability (Phonons)

6. Tools for producing comprehensive reports, convex hulls, band structures,
   density of states, etc

7. Datamining tools to extract knowledge from the structures found, identify
   patterns in the data and identify suitable candidates for technological
   applications

8. A web interface


This is a new project and many classes and methods are refactored frequently.
This is and will be a work in progress. We hope to stabilize the most critical
classes for the release 1.

This code is open-source. We also welcome extra hands to improve this library
with your own contributions. At present only one developer has being in charge
of the project. More hands and eyes are very welcomed.

Contributors
------------

1. Prof. Aldo H. Romero [West Virginia University] (Project Director)

1. Guillermo Avendaño-Franco [West Virginia University] (Basic Infraestructure)

1. Irais Valencia Jaime [West Virginia University] (Simulation and testing)

1. Sobhit Singh [West Virginia University] (Data-mining)

1. Francisco Muñoz (Universidad de Chile)