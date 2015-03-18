Database (pychemia.db)
======================

PyChemia stores the structures in a MongoDB database. The database
contains several collections, such as

* pychemia_entries: The main collection, contains all the structures,
properties and status info

* fingeprints: Collection with using the same IDs as pychemia_entries
storing the fingerprints of the corresponding structures

* generations: Collection with the same IDs as pychemia_entries to
store on which generation each structure belongs.

* generation_changes: Collection storing the changes introduced from one
generation to the next one.

* population_info: Stores the basic parameters entered to create the population

* searcher_info: Store parameters about the searcher use to populate the database

.. automodule:: pychemia.db
   :members:
