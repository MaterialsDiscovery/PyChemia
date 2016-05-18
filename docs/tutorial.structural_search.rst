PyChemia Software Framework
===========================

PyChemia is framework for materials discovery, more than just compute
DFT calculations for structures already present in databases such as ICSD,
PyChemia searches for new structures that could never have been synthesized
or reported before.

To achieve its goal, PyChemia relies on a number of methods of structural
search such as minima hoping method, genetic algorithms and other population
based algorithms.

As a software framework, PyChemia is structurated around five axis.
They are:

1. Structural Search Methods

2. Storage and Databases

3. Data Mining or Knowledge Discovery in Data

4. Execution and Analysis

5. Reporting and Visualization


We will develop those five axis and the relations between them

Structural Search Methods
-------------------------

One of the distinctive characteristics of PyChemia is its ability to search for new structures.
Ab-initio calculations of structures present of databases such as ICSD becomes limited to the extension of the original
databases.
Even if such effort is indeed a big challenge with large databases, those databases represent a small portion of the
structures that could be created.
The ability to predict new structures and target those findings to specific applications is a task of technological
relevance.

PyChemia was created with an strong focus on structural prediction.
PyChemia implements several methods, from minima hoping method to metaheuristics


Minima Hoping Method
~~~~~~~~~~~~~~~~~~~~


Metaheuristics
~~~~~~~~~~~~~~


Storage and Databases
---------------------

If PyChemia were a software to compute ab-initio calculations from structures taken from another database or from
very predictive set of prototypes only one database will be enough.
However, as we mention before PyChemia is focused on structural search.

The kind of algorithms that we describe above typically explore thousands of different structures before selecting a
reduced subset of thermodynamically stable or metastable ones.

Also, the search for new structures must be guided by specific applications of interest, batteries, thermoelectrics,
superconductors, etc. Different applications are associated to different properties in the electronic structure.

Those two elements, the dynamic nature of structural search methods and the flexibility in the data that we would
like to store is the reason why we are not using a single database and why we departed from traditional SQL schemas.

PyChemia was designed to create new databases for each structural search. We still keep one large database with a
curated selection of structures but the ability to create small and flexible databases is central for the success
of PyChemia.

PyChemia relies on MongoDB. MongoDB is an open-source document database, and the leading NoSQL database.
We take advantage of dynamic schemas that offer simplicity and power and are in harmony with the own principles
of scientific research. Different structures are intended for different applications and different applications
requieres the calculation of different physical properties.

Data Mining or Knowledge Discovery in Data
------------------------------------------

Data Mining or more correctly Knowledge Discovery in Data is the computational process of discovering patterns in
large data sets involving methods at the intersection of artificial intelligence, machine learning, statistics,
and database systems.

Structural Search algorithms has the ability to gererate in the long term more structures that those that could
be synthesized in reasonable time.
As the database becomes large traditional techniques based on estimate good candidates based on just a few properties
will be replaced by automatize algorithms that search patterns in structures and predict were new stochiometries
deserve exploration.

Bayesian networks and Gaussian Processes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Reporting and Visualization
---------------------------

Doing electronic structure calculations for thousands of materials is challenging but there is not point in doing that
without the ability to communicate to experimentalist structures that worth to be synthesized.
PyChemia uses Django as a web-frontend to display and explore the computations not only on the global database but also
on the specific database product of a structural search run.

Execution and Analysis
----------------------

There are two levels of execution in structural search runs. Ab-initio calculations and structural analysis.
PyChemia provides concurrent execution and basic job management for running multiple ab-initio calculations
on clusters or non-queue systems

PyChemia relies on state-of-the-art ab-initio software packages to compute electronic structure and their properties.
In particular PyChemia has support for VASP, ABINIT, Octopus, Fireball and DFTB+

Different levels of theory allows for a better compromise between accuracy and computational cost.

Structural analysis is defined as all kinds of procedures that relies only on structural data, for example
structural fingerprints, topological analysis, hardness calculations, comparators and symmetry calculations.
PyChemia provides multi-threaded execution of structural analysis routines


