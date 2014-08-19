from pymongo import MongoClient
import collections

import pychemia


__author__ = 'Guillermo Avendano Franco'


class PyChemiaDB():

    def __init__(self, name='pychemiadb', host='localhost', port=27017):

        client = MongoClient(host, port)
        db = client[name]

        self.db = db
        self.structures_col = self.db.structures
        self.properties_col = self.db.properties

    def insert(self, structure, properties):

        struct_id = self.structures_col.insert(structure)
        properties['struct_id'] = struct_id
        prop_id = self.properties_col.insert(properties)
        return struct_id, prop_id

    def get_iterator(self):
        cursor = self.db.structures.find()
        return Iterator(self.db, cursor)


class Iterator(collections.Iterable):

    def __init__(self, db, cursor):
        self.db = db
        self.cursor = cursor

    def __iter__(self):
        return self

    def __next__(self):
        structure_entry = self.cursor.next()
        structure = pychemia.Structure().fromdict(structure_entry)
        properties = self.db.properties.find_one({'struct_id': structure_entry['_id']})
        return structure, structure_entry, properties

    next = __next__
