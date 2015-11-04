"""
Routines related to Metadata info and Repositories
"""

try:
    import pymongo

    if pymongo.version_tuple[0] < 3:
        print 'Pymongo version is too old for PyChemia, install pymongo 3+'
        USE_MONGO = False
    else:
        USE_MONGO = True
except ImportError:
    pymongo = None
    print 'Could no import pymongo, mongo database functionality disabled'
    USE_MONGO = False

try:
    import gridfs

    USE_GRIDFS = True
except ImportError:
    gridfs = None
    print 'Could no import pymongo, mongo database functionality disabled'
    USE_GRIDFS = False

if USE_MONGO:
    from _db import PyChemiaDB, get_database, object_id, create_user, create_database

    if USE_GRIDFS:
        from _queue import PyChemiaQueue

# __all__ = filter(lambda s: not s.startswith('_'), dir())
