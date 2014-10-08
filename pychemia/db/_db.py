__author__ = 'Guillermo Avendano Franco'

from pymongo import MongoClient
import collections
import pychemia


class PyChemiaDB():

    def __init__(self, name='pychemiadb', host='localhost', port=27017, master=False):

        self.name = name
        self._is_master = master
        self._client = MongoClient(host, port)
        self.db = self._client[name]
        self.entries = self.db.pychemia_entries

    def insert(self, structure, properties=None):
        """
        Insert a pychemia structure instance and properties
        into the database
        :param structure: (pychemia.Structure) An instance of Pychemia's Structure
        :param properties: (dict) Dictionary of properties
        :return:
        """
        entry_dict = structure.todict()
        if properties is not None:
            entry_dict['properties'] = properties
        else:
            entry_dict['properties'] = {}
        entry_id = self.entries.insert(entry_dict)
        return entry_id,

    def get_iterator(self):
        cursor = self.db.structures.find()
        return Iterator(self.db, cursor)

    def clean(self):
        self._client.drop_database(self.name)
        self.db = self._client[self.name]

    def is_master(self):
        return self._is_master

