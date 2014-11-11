"""
Routines related to Metadata info and Repositories
"""

from _repo import StructureEntry, ExecutionRepository, PropertiesEntry
try:
    from _db import PyChemiaDB
    USE_MONGO = True
except ImportError:
    print 'PYMONGO not found, mongo database functionality disabled'
    USE_MONGO = False


#__all__ = filter(lambda s: not s.startswith('_'), dir())

