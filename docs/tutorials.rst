Basic Tutorial on PyChemia
==========================

The core of PyChemia as a library is the Structure class.
The Structure class defines atomic structures that could be periodic or not
and with periodicity defined in 1, 2 or 3 directions.
This tutorial will show the basic functionality of the PyChemia Structure 
class.

This document offers and introduction to the use of the library.
This information is useful for both users and developers to have get
a landscape about how the code is structured.

PyChemia uses strognly the Object oriented capabilities of Python.
This code is not a set of scripts with functions. Almost everything
is encapsulated in a class.

The code was build on a first instance for crystal structure prediction,
but the Structure object is general enough to support atomic structures
with no periodicity (moelcules and clusters), periodicity in one direction
(wires and tubes), periodicity in two directions (thin films and graphene like
structures) and fully periodic structures.

There are two basic modules that are the backbone of the entire library
pychemia.geometry and pychemia.repositories

1. :mod:`pychemia.geometry`: Contains the class Structure that stores all the
   geometric information about a crystal or molecule and provides methods to
   manipulate its contents. This class must be highly flexible because the
   purpose of structural search is precisely change structures to find new
   ones.

   The Structure object contains a description of atoms and cell parameters.
   Its main variables are:

   Structure.symbols: The set of symbols in the cell
                      Example [ 'Na', 'Cl']

   Structure.positions: Dimensional cartesian coordinates

   Structure.reduced: non-dimensional cell reduced coordinates
                      Example [[0,0],[0, 0.5]]

   Structure.natom: Number of atoms in the cell

   Structure.cell: The lattice vectors


2. :mod:`pychemia.repositories`: This module contains four important classes
   to manipulate collections of structures. StructureRepository defines the
   entries in a directory collecting many structures. This class keep an eye
   on each structure (StructureEntry) that is added or removed from the
   repository. This object is also able to merge another StructureRepository
   and search for eventual duplicates.
   A ExecutionRepository for the other hand stores Executions this is a heavy
   directory where every single calculation is stored. Calculations should be
   deleted only if the structure for which this calculation was done is also
   removed from the main repository.


Basic Workflow for PyChemia
---------------------------

1. We receive a set of CIF files to add in the database, with a set of
   different compositions

2. The code read the CIF files, create Structure objects and feed the repositories
   with them

3. The selected structural predictor takes those structures or create new ones based
   on the method of choice.
   The method returns new Structures those are stored in a separate repository, they
   are screened by similarity with those already present.

4. The structures enter in a set of 'tasks' for thermal stability and dynamic stability.

5. Survival structures are stored in the final repository and final characterizations
   can follow on those Structures

