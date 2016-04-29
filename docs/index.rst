.. PyChemia documentation master file

PyChemia
========

PyChemia is a python library for automatize atomistic simulations.
PyChemia was created with a particular emphasis on materials
discovery an analysis using first principles approaches.
The idea behind the library is provide basic building blocks that
allow reasercher create complex tasks with minimal scripting.

You can consider PyChemia as library build in a set of layers of
functionality.

Tier 0:

The core module of PyChemia is made around three classes: Structure,
Composition and Lattice. Almost all modules work on top of those
three classes in some way.

Tier 1:

There are a set of modules that offer the functionality required for
complex tasks of atomistic simulations. The classes are "db" for
database interaction, "symm" for symmetry analysis, "code" for
interaction with external atomistic calculators and "analysis"
that offer classes for geometric analysis of structures.

Tier 2:

It is common for complex tasks on atomistic simuations to deal with
a population of candidades, one example is "global search".
The module population implement the functionality to deal with a set
of structures stored on the database and comunicate when a given entry
is evaluated or not.

Tier 3:

PyChemia original conception was around structural search, a task that
involves using a global searcher to find competitive structures by
exploring the potential energy surface. This original functionality
has been largely generalize to include global search beyond structural
search. 


Contents:

.. toctree::

   intro
   core
   tutorial
   db
   structural_search_general
   structural_search_pychemia

   
Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

