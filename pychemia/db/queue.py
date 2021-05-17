import os
import pymongo
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
                uri += ':' + str(passwd)
            uri += '@'
        uri += host + ':' + str(port)
        print('URI:', uri)
        if user is not None:
            uri += '/' + name
        if replicaset is not None:
            self._client = pymongo.MongoClient(uri, ssl=ssl, replicaset=replicaset)
        else:
            self._client = pymongo.MongoClient(host=host, port=port, ssl=ssl, tlsAllowInvalidCertificates=True)
        for i in ['version']:
            print('%20s : %s' % (i, self._client.server_info()[i]))
        self.db = self._client[name]
        if user is not None and self.db.authenticate(user, passwd):
            print('Authentication successful')

        self.set_minimal_schema()
        self.fs = gridfs.GridFS(self.db)

    def add_file(self, entry_id, location, filepath):

        assert (os.path.isfile(filepath))
        hashcode = hashfile(filepath)
        rf = open(filepath, 'rb')
        filename = os.path.basename(filepath)
        length = os.path.getsize(filepath)

        existing = self.db.fs.files.find_one({'hash': hashcode, 'length': length, 'filename': filename})

        if existing is None:

            file_id = self.fs.put(rf, filename=os.path.basename(filename), hash=hashcode)
            print('New file ', file_id)
            self.db.pychemia_entries.update_one({'_id': entry_id}, {'$addToSet': {location + '.files': {'file_id': file_id,
                                                                                                    'name': filename,
                                                                                                    'hash': hashcode}}})
        else:
            file_id = existing['_id']
            print('File already present ', file_id)
            self.db.pychemia_entries.update_one({'_id': entry_id}, {'$addToSet': {location + '.files': {'file_id': file_id,
                                                                                                    'name': filename,
                                                                                                    'hash': hashcode}}})

    def add_input_file(self, entry_id, filename):
        self.add_file(entry_id, 'input', filename)

    def set_minimal_schema(self):
        for entry in self.db.pychemia_entries.find({'meta': None}, {'_id': 1}):
            print('Missing field "meta" on', entry['_id'])
            self.db.pychemia_entries.update_one({'_id': entry['_id']}, {'$set': {'meta': {}}})
        for entry in self.db.pychemia_entries.find({'input': None}, {'_id': 1}):
            print('Missing field "input" on', entry['_id'])
            self.db.pychemia_entries.update_one({'_id': entry['_id']}, {'$set': {'input': {}}})
        for entry in self.db.pychemia_entries.find({'output': None}, {'_id': 1}):
            print('Missing field "output" on', entry['_id'])
            self.db.pychemia_entries.update_one({'_id': entry['_id']}, {'$set': {'output': {}}})
        for entry in self.db.pychemia_entries.find({'job': None}, {'_id': 1}):
            print('Missing field "job" on', entry['_id'])
            self.db.pychemia_entries.update_one({'_id': entry['_id']}, {'$set': {'job': {}}})

    def set_structure(self, entry_id, location, structure):
        assert (location in ['input', 'output'])
        self.db.pychemia_entries.update_one({'_id': entry_id}, {'$set': {location + '.structure': structure.to_dict}})

    def set_input(self, entry_id, code, inputvar):

        for i in inputvar:
            if i.startswith('$'):
                value = inputvar.pop(i)
                inputvar[i[1:]] = value
        print('INPUTVAR:\n %s %s' % (inputvar, dict(inputvar)))
        self.db.pychemia_entries.update_one({'_id': entry_id}, {'$set': {'input.variables': dict(inputvar),
                                                                     'input.code': code.lower()}})

    def new_entry(self, structure=None, variables=None, code=None, files=None, priority=0, dbname=None, db_id=None):

        if variables is not None and code is None:
            raise ValueError("Input variables requiere code name")
        if variables is None and code is not None:
            raise ValueError("Input variables require code name")

        entry = {'input': {}, 'output': {}, 'job': {}, 'meta': {'submitted': False,
                                                                'priority': priority,
                                                                'finished': False,
                                                                'deployed': False,
                                                                'dbname': dbname,
                                                                'db_id': db_id}}
        result = self.db.pychemia_entries.insert_one(entry)
        entry_id = result.inserted_id

        # Commented for compatibility with mongo 2.4
        # self.db.pychemia_entries.update_one({'_id': entry_id}, {'$currentDate': {'meta.CreationDate': True}})

        if structure is not None:
            self.set_input_structure(entry_id, structure)
        if variables is not None and code is not None:
            self.set_input(entry_id, code=code, inputvar=variables)
        if files is not None:
            for ifile in files:
                self.add_input_file(entry_id, filename=ifile)
        return entry_id

    def set_job_settings(self, entry_id, nparal=None, queue=None, nhours=None, mail=None, task_name=None,
                         task_settings=None, task_kind=None):
        if nparal is not None:
            self.db.pychemia_entries.update_one({'_id': entry_id}, {'$set': {'job.nparal': nparal}})
        if queue is not None:
            self.db.pychemia_entries.update_one({'_id': entry_id}, {'$set': {'job.queue': queue}})
        if mail is not None:
            self.db.pychemia_entries.update_one({'_id': entry_id}, {'$set': {'job.mail': mail}})
        if nhours is not None:
            self.db.pychemia_entries.update_one({'_id': entry_id}, {'$set': {'job.nhours': nhours}})
        if task_name is not None:
            self.db.pychemia_entries.update_one({'_id': entry_id}, {'$set': {'job.task_name': task_name}})
        if task_kind is not None:
            self.db.pychemia_entries.update_one({'_id': entry_id}, {'$set': {'job.task_kind': task_name}})
        if task_settings is not None:
            self.db.pychemia_entries.update_one({'_id': entry_id}, {'$set': {'job.task_settings': task_settings}})

    def set_input_structure(self, entry_id, structure):
        return self.set_structure(entry_id, 'input', structure)

    def set_output_structure(self, entry_id, structure):
        return self.set_structure(entry_id, 'output', structure)

    def get_input_variables(self, entry_id):
        entry = self.db.pychemia_entries.find_one({'_id': entry_id}, {'input.variables': 1})
        if 'variables' in entry['input']:
            return entry['input']['variables']
        else:
            return None

    def get_structure(self, entry_id, location):
        """
        Return the structure in the entry with id 'entry_id'

        :rtype : Structure
        """
        assert (location in ['input', 'output'])
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
