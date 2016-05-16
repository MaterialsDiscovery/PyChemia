"""
PyChemia relies of MongoDB to store structures and the properties computed by atomistic codes.
This module creates and manipulates Mongo databases.
There are two kinds of databases, PyChemiaDB is a kind of database to store structure and properties
PyChemiaQueue is a repository of calculations.
"""

from pychemia import HAS_PYMONGO, HAS_GRIDFS

if HAS_PYMONGO:
    from .db import PyChemiaDB, get_database, object_id, create_user, create_database

    if HAS_GRIDFS:
        from .queue import PyChemiaQueue

# __all__ = filter(lambda s: not s.startswith('_'), dir())
