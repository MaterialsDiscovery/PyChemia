"""
Routines related to Metadata info and Repositories
"""

from _repo import StructureEntry, ExecutionRepository, PropertiesEntry
try:
    from _db import PyChemiaDB, get_database, object_id
    USE_MONGO = True
except ImportError:
    print 'Could no import pymongo, mongo database functionality disabled'
    USE_MONGO = False


# __all__ = filter(lambda s: not s.startswith('_'), dir())

