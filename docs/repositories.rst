Repositories (pychemia.repositories)
====================================

PyChemia stores the structures in a directory structure called a
repository, each structure is identified with a UUID string
and the corresponding directory contains at least three files

* structure.json : The serialized content of the structure object
* properties.json : A dictionary with the properties computed for
  that structure
* metadata.json : Information about tags, parents and children of the
  structure

.. automodule:: pychemia.db
   :members:


Structure Entry
---------------

.. autoclass:: pychemia.db.StructureEntry
   :members:

Properties Entry
----------------

.. autoclass:: pychemia.db.PropertiesEntry
   :members:

Structure Repository
--------------------

.. autoclass:: pychemia.db.StructureRepository
   :members:

Execution Repository
--------------------

.. autoclass:: pychemia.db.ExecutionRepository
   :members:

