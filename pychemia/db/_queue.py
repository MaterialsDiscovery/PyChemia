__author__ = 'Guillermo Avendano-Franco'

import os
from pymongo import MongoClient
from pychemia import Structure
import gridfs
from pychemia.utils.computing import hashfile

class PyChemiaQueue:

    def __init__(self, name='Queue', host='localhost', port=27017, user=None, passwd=None, ssl=False, replicaset=None):
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
        if replicaset is not None:
            self._client = MongoClient(uri, ssl=ssl, replicaset=replicaset)
        else:
            self._client = MongoClient(uri, ssl=ssl)
        self.db = self._client[name]
        self.set_minimal_schema()
        self.fs = gridfs.GridFS(self.db)

    def add_file(self, entry_id, location, filepath):

        assert(os.path.isfile(filepath))
        hashcode = hashfile(filepath)
        rf = open(filepath, 'rb')
        filename = os.path.basename(filepath)
        length = os.path.getsize(filepath)

        existing = self.db.fs.files.find_one({'hash': hashcode, 'length': length, 'filename': filename})

        if existing is None:

            file_id = self.fs.put(rf, filename=os.path.basename(filename), hash=hashcode)
            print 'New file ', file_id
            self.db.pychemia_entries.update({'_id': entry_id}, {'$addToSet': {location+'.files': {'file_id': file_id,
                                                                                                  'name': filename,
                                                                                                  'hash': hashcode}}})
        else:
            file_id = existing['_id']
            print 'File already present ', file_id
            self.db.pychemia_entries.update({'_id': entry_id}, {'$addToSet': {location+'.files': {'file_id': file_id,
                                                                                                  'name': filename,
                                                                                                  'hash': hashcode}}})


    def add_input_file(self, entry_id, filename):
        self.add_file(entry_id, 'input', filename)

    def set_minimal_schema(self):
        for entry_id in self.db.pychemia_entries.find({'meta': None}, {'meta': 1}):
            print entry_id
            self.db.pychemia_entries.update({'_id': entry_id['_id']}, {'$set': {'meta': {}}})
        for entry_id in self.db.pychemia_entries.find({'input': None}):
            print entry_id
            self.db.pychemia_entries.update({'_id': entry_id['_id']}, {'$set': {'input': {}}})
        for entry_id in self.db.pychemia_entries.find({'output': None}):
            print entry_id
            self.db.pychemia_entries.update({'_id': entry_id['_id']}, {'$set': {'output': {}}})

    def set_structure(self, entry_id, location, structure):
        assert(location in ['input', 'output'])
        self.db.pychemia_entries.update({'_id': entry_id}, {'$set': {location+'.structure': structure.to_dict}})

    def set_input(self, entry_id, code, input):
        self.db.pychemia_entries.update({'_id': entry_id}, {'$set': {'input.variables': input.variables,
                                                                     'input.code': code.lower()}})

    def new_entry(self, structure=None, variables=None, code=None, files=None):

        if variables is not None and code is None:
            raise ValueError("Input variables requiere code name")
        if variables is None and code is not None:
            raise ValueError("Input variables requiere code name")

        entry = {'input': {}, 'output': {}, 'meta': {}}
        entry_id = self.db.pychemia_entries.insert(entry)

        if structure is not None:
            self.set_input_structure(entry_id, structure)
        if variables is not None and code is not None:
            self.set_input(entry_id, code=code, input=variables)
        if files is not None:
            for ifile in files:
                self.add_input_file(entry_id, filename=ifile)
        return entry_id

    def set_input_structure(self, entry_id, structure):
        return self.set_structure(entry_id, 'input', structure)

    def set_output_structure(self, entry_id, structure):
        return self.set_structure(entry_id, 'output', structure)

    def get_input_variables(self, entry_id):
        entry = self.db.pychemia_entries.find_one({'_id': entry_id}, {'input.variables': 1})
        return entry['input']['variables']

    def get_structure(self, entry_id, location):
        """
        Return the structure in the entry with id 'entry_id'

        :rtype : Structure
        """
        assert(location in ['input', 'output'])
        entry = self.db.pychemia_entries.find_one({'_id': entry_id}, {location: 1})
        return Structure.from_dict(entry[location]['structure'])

    def get_input_structure(self, entry_id):
        return self.get_structure(entry_id, 'input')

    def get_output_structure(self, entry_id):
        return self.get_structure(entry_id, 'output')

    def write_input_files(self, entry_id, destination=None):

        if destination is None:
            dest = '.'
        elif os.path.isfile(destination):
            dest = os.path.dirname(os.path.abspath(destination))
        elif not os.path.exists(destination):
            os.mkdir(destination)
            dest = destination
        elif os.path.isdir(destination):
            dest = destination
        else:
            raise ValueError('Destination not valid')

        entry = self.db.pychemia_entries.find_one({'_id': entry_id}, {'input.files': 1})
        for ifile in entry['input']['files']:
            rf = self.fs.get(ifile['file_id'])
            wf = open(dest + os.sep + ifile['name'], 'wb')
            wf.write(rf.read())
            wf.close()
