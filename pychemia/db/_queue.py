__author__ = 'Guillermo Avendano-Franco'

import os
import hashlib
from pymongo import MongoClient
from pychemia import Structure
import gridfs


class PyChemiaQueue:

    def __init__(self, name='pychemiadb', host='localhost', port=27017, user=None, passwd=None, ssl=False,
                 replicaset=None):
        """
        Creates a MongoDB client to 'host' with 'port' and connect it to the database 'name'.
        Authentication can be used with 'user' and 'password'

        :param name: (str) The name of the database
        :param host: (str) The host as name or IP
        :param port: (int) The number of port to connect with the server (Default is 27017)
        :param user: (str) The user with read or write permissions to the database
        :param passwd: (str/int) Password to authenticate the user into the server

        :return:
        """
        self.db_settings = {'name': name, 'host': host, 'port': port, 'user': user, 'passwd': passwd, 'ssl': ssl}
        self.name = name
        uri = 'mongodb://'
        if user is not None:
            uri += user
            if passwd is not None:
                uri += ':'+str(passwd)
            uri += '@'
        uri += host + ':' + str(port)
        if user is not None:
            uri += '/' + name
        self._client = MongoClient(uri, ssl=ssl, replicaset=replicaset)
        self.db = self._client[name]
        self.set_minimal_schema()
        self.fs = gridfs.GridFS(self.db)

    def add_file(self, entry_id, location, filename):

        assert(os.path.isfile(filename))
        hashcode = hashfile(filename)
        rf = open(filename, 'rb')
        file_id = self.fs.put(rf, filename=os.path.basename(filename), hash=hashcode)
        self.db.pychemia_entries.update({'_id': entry_id}, {'$addToSet': {location+'files': file_id}})

    def set_minimal_schema(self):
        for entry_id in self.db.pychemia_entries.find({'status': None}, {'status': 1}):
            print entry_id
            self.db.pychemia_entries.update({'_id': entry_id['_id']}, {'$set': {'status': {}}})
        for entry_id in self.db.pychemia_entries.find({'input': None}):
            print entry_id
            self.db.pychemia_entries.update({'_id': entry_id['_id']}, {'$set': {'input': {}}})
        for entry_id in self.db.pychemia_entries.find({'output': None}):
            print entry_id
            self.db.pychemia_entries.update({'_id': entry_id['_id']}, {'$set': {'output': {}}})

    def set_structure(self, location, structure, entry_id):
        assert(location in ['input', 'output'])
        self.db.pychemia_entries.update({'_id': entry_id}, {'$set': {location+'.structure': structure.to_dict}})

    def new_entry(self):
        entry = {'input': {}, 'output': {}, 'status': {}}
        entry_id = self.db.pychemia_entries.insert(entry)
        return entry_id

    def set_input_structure(self, structure, entry_id):
        return self.set_structure('input', structure, entry_id)

    def set_output_structure(self, structure, entry_id):
        return self.set_structure('output', structure, entry_id)

    def get_structure(self, location, entry_id):
        """
        Return the structure in the entry with id 'entry_id'

        :rtype : Structure
        """
        assert(location in ['input', 'output'])
        entry_id = object_id(entry_id)
        entry = self.db.pychemia_entries.find_one({'_id': entry_id}, {location: 1})
        return Structure.from_dict(entry[location]['structure'])

    def get_input_structure(self, entry_id):
        return self.get_structure('input', entry_id)

    def get_output_structure(self, entry_id):
        return self.get_structure('output', entry_id)


def hashfile(filename):
    blocksize = 65536
    hasher = hashlib.md5()
    with open(filename, 'rb') as afile:
        buf = afile.read(blocksize)
        while len(buf) > 0:
            hasher.update(buf)
            buf = afile.read(blocksize)
    return hasher.hexdigest()