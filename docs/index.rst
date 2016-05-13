.. PyChemia documentation master file

PyChemia
========

PyChemia is a python library for automatize atomistic simulations.
PyChemia is build around a core module with two classes and a set of
some other modules that offers a variety of operations in order to
perform more complex operations.

The core is made of two classes 'Structure' and 'Composition'.
Other modules are:

   * analysis: Uses pure structure information for changing structures,
   	matching atoms between two structures, create slabs, surfaces and some
   	other geometrical operations on structures.

   * code: This module deals with creating inputs and reading outputs from several atomistic
   	simulation codes. Right now, we support ABINIT, DFTB+, Fireball, an internal
   	calculator for LennardJones Clusters, Octopus and VASP.

   * core: The core of PyChemia are two classes that are imported at the root level of the library:
   	Structure and Composition. For PyChemia, Structure is a set of sites with one or more atoms located
   	on each site with a probability associated to them. The Structure could be finite or periodic
   	in one or more directions. In the case of a crystal the Structure will have also a lattice.
   	Composition is a set o atoms of a define species. For crystals is the set of atoms on each unit cell.
   	No geometry is store on a Composition object, and the order of atoms is irrelevant.

   * crystal: For the case of periodic structures in three directions, a set of modules is created
   	for three basic properties of a crystal. The class 'KPoints' store the description a a k-point
   	mesh, path, or direct list of points on the reciprocal lattice. The class 'Lattice' store and
   	manipulate the cell vectors of a crystal. The third class is 'CrystalSymmetry' for computing
   	the Space Group, getting structures on the Bravais cell and finding the primitive cell.

   * db: PyChemia uses MongoDB as Database for storing Structures and the properties computed
   	by atomistic simulation codes. There are two classes defined, PyChemiaDB to store structures
   	and properties and PyChemiaQueue to store a massive set of calculations and the status of those
   	calculations.

   * dm: This is a module in development, it contains classes for DataMining PyChemia Databases
   	and global searches. Right now it contains a class for Network analysis, but in the future
   	will provide interfaces with some other libraries for machine learning.

   * evaluator: There are two circumstances where a atomistic simulation can be perform. If you
   	have a machine without a queue system PyChemia will provide a very simplified queue for computing
   	concurrent calculations under the constrains of your number of cores. If you use a queue
   	system, the evaluator will setup the batch scripts, submit the jobs and monitor the status
   	of those jobs, finally when the job is done, it will collect the final data and update the
   	the databases.

   * external: There are some other packages for which PyChemia worth interact, ASE (Atomistic
   	Simulation Environment) is a python package supporting a remarkable  number of calculators.
   	Findsym is a close software to compute spacegroups and computing the CIF file of a given
   	structure. 'pymatgen' is python code behind 'Materials Project' and a outstanding piece of
   	very well written code. Originally implemented for VASP only now includes support for ABINIT.

   * io: There are a pletora of formats for describing atomistic structures. We support three
   	basic file-formats, ascii, CIF and XYZ. This module provides the classes for reading and
   	writing them.

   * md: Molecular Dynamics is an important operation for atomistic simulations. This module
   	provide a 'in-house' calculator for MD simulations, using the forces computed from static
   	calculations from an external atomistic code.

   * plotting: PyChemia includes a set of classes and routines for interfacing with some other
   	libraries and external software for data visualization, 3D plotting and graphic representation.
   	pyprocar plots band-structures, LatticePlot and StructurePlot uses mayavi for visualizing
   	atomic configurations and lattices. We have also interfaces for Povray and a developing interface
   	to D3.js.

   * population: A population is basically a collection of candidates collected 
	somehow by a global searcher or from a DB query. PyChemia includes 
	classes for populations of Lennard-Jones clusters, populations of 
	vectors on a N-dimensional space.

   * runner: Controls the executions in cluster by creating batch scripts 
	and checking the status of the queue. Right now only supports Torque

   * searcher: Global search operations, the prototypical case is structural 
	search but any pychemia population could be use for searching. Several
     	metaheuristic algorithms have been implemented

   * utils: Several small routines, like a periodic table, mathematical 
	operations, conversions and serializers. 

   * web: On development, creation of a web interface for looking into the 
	pychemia databases and controling executions.


Contents:

.. toctree::

   intro
   core
   analysis
   code
   crystal
   db
   dm
   evaluator
   external
   io
   md
   population
   runner
   searcher
   utils
   visual
   web
   tutorials
   structural_search_pychemia
   structural_search_general

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
