pychemia.db
===========

One core capability of PyChemia is global optimization of atomic
structures such as clusters and crystals.
Doing global optimization requires compute forces and energies of
hundreds or even thousands of structures and all these data must
be store and processed efficiently.

PyChemia stores the structures and their properties in a MongoDB database. 
The database contains several collections, such as

* pychemia_entries: The main collection, contains all the structures,
                    properties and status info

* fingeprints: Collection with using the same IDs as pychemia_entries
               storing the fingerprints of the corresponding structures

* generations: Collection with the same IDs as pychemia_entries to
               store on which generation each structure belongs.

* generation_changes: Collection storing the changes introduced from one
                      generation to the next one.

* population_info: Stores the basic parameters entered to create the population

* searcher_info: Store parameters about the searcher use to populate 
                 the database

* lineage: Store the heritage of structures produced by the global searcher

To do a more efficient executiion of calculations using several computer clusters
and dedicated machines, PyChemia provides a central database of calculations
that can be used to feed clusters with execution jobs. 
The PyChemiaQueue is used to fill the role of a meta-queue for structures and 
jobs that need to be computed. 

.. automodule:: pychemia.db
:members:
