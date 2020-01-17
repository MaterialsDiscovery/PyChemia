import itertools
import json
import socket
from multiprocessing import Pool

import numpy as np
import pymongo
from pymongo.errors import ServerSelectionTimeoutError, OperationFailure
from bson.objectid import ObjectId

from pychemia import Structure, HAS_PYMONGO
from pychemia.utils.periodic import atomic_symbols


class PyChemiaDB:
    def __init__(self, name='pychemiadb', host='localhost', port=27017, user=None, passwd=None, ssl=False,
                 replicaset=None):
        """
        Creates a MongoDB client to 'host' with 'port' and connect it to the database 'name'.
        Authentication can be used with 'user' and 'password'

        :param name: (str) The name of the database
        :param host: (str) The host as name or IP
        :param port: (int) The number of port to connect with the server (Default is 27017)
        :param user: (str) The user with read or write permissions to the database
        :param passwd: (str,int) Password to authenticate the user into the server

        """
        self.db_settings = {'name': name,
                            'host': host,
                            'port': port,
                            'user': user,
                            'passwd': passwd,
                            'ssl': ssl,
                            'replicaset': replicaset}
        self.name = name
        maxSevSelDelay = 2
        uri = 'mongodb://'
        if user is not None:
            uri += user
            if passwd is not None:
                uri += ':' + str(passwd)
            uri += '@'
        uri += host + ':' + str(port)
        if user is not None:
            uri += '/' + name
        try:
            if pymongo.version_tuple[0] == 2:
                self._client = pymongo.MongoClient(uri, ssl=ssl, replicaSet=replicaset,
                                                   serverSelectionTimeoutMS=maxSevSelDelay)
            elif pymongo.version_tuple[0] == 3:
                self._client = pymongo.MongoClient(uri, ssl=ssl,
                                                   ssl_cert_reqs=pymongo.ssl_support.ssl.CERT_NONE,
                                                   replicaSet=replicaset, serverSelectionTimeoutMS=maxSevSelDelay)
            else:
                raise ValueError('Wrong version of pymongo')

            try:
                self._client.server_info()
            except OperationFailure:
                raise RuntimeError("ERROR: Database '%s' on '%s' cannot be accessed, either the database does not "
                                   "exist or you need the right credentials to get access" % (name, host))

        except ServerSelectionTimeoutError:
            raise RuntimeError("ERROR: No connexion could be established to server: %s" % uri)

        self.db = self._client[name]
        self.entries = self.db.pychemia_entries
        self.set_minimal_schema()

    def __str__(self):
        ret = ' Database Name:       %s\n' % self.name
        ret += ' Host:                %s\n' % self.db_settings['host']
        ret += ' Port:                %s\n' % self.db_settings['port']
        ret += ' User:                %s\n' % self.db_settings['user']
        ret += ' SSL:                 %s\n' % self.db_settings['ssl']
        return ret

    def save_json(self, filename='db_settings.json'):
        wf = open(filename, 'w')
        json.dump(self.db_settings, wf, sort_keys=True, indent=4, separators=(',', ': '))
        wf.close()

    def insert(self, structure, properties=None, status=None, entry_id=None):
        """
        Insert a pychemia structure instance and properties
        into the database
        :param entry_id: Mongo ID for the entry
        :param structure: (pychemia.Structure) An instance of Pychemia's Structure
        :param properties: (dict) Dictionary of properties
        :param status: (dict) Dictionary of status
        :return:
        """
        if status is None:
            status = {}
        if properties is None:
            properties = {}
        entry = {'structure': structure.to_dict, 'properties': properties, 'status': status}
        if entry_id is not None:
            entry['_id'] = entry_id
        result = self.entries.insert_one(entry)
        return result.inserted_id

    def clean(self):
        self._client.drop_database(self.name)
        self.db = self._client[self.name]

    def update(self, entry_id, structure=None, properties=None, status=None):
        """
        Update the fields 'structure', 'properties' or 'status' for a given identifier 'entry_id'

        :param entry_id: (ObjectID, str)
        :param structure: (pychemia.Structure) Structure to update
        :param properties: (dict) Dictionary of properties to update
        :param status: (dict) Status dictionary

        :return: The identifier for the entry that was updated

        :rtype : ObjectId

        """

        assert (self.entries.find_one({'_id': entry_id}) is not None)
        entry = self.entries.find_one({'_id': entry_id})
        if structure is not None:
            if isinstance(structure, Structure):
                entry['structure'] = structure.to_dict
            elif isinstance(structure, dict):
                entry['structure'] = structure
            else:
                print('ERROR: Could not process the structure')
                print(type(structure))
        if properties is not None:
            entry['properties'] = properties
        if status is not None:
            entry['status'] = status

        self.entries.replace_one({'_id': entry_id}, entry)
        return entry_id

    def find_AnBm(self, specie_a=None, specie_b=None, n=1, m=1):
        """
        Search for structures with a composition expressed as AnBm
        where one and only one between A or B is fixed and the numbers
        amounts n and m are both fixed

        :param specie_a: (str) atom symbol for the first specie
        :param specie_b: (str) atom symbol for the second specie
        :param n: number of atoms for specie 'a'
        :param m: number of atoms for specie 'b'
        :return:
        :return: (list) List of ids for all the structures that fulfill
                 the conditions
        """
        if specie_a is None and specie_b is None:
            raise ValueError("Enter a specie for A or B")
        elif specie_a is not None and specie_b is not None:
            raise ValueError("Only enter A or B, not both")
        elif specie_a is not None:
            atom_fixed = specie_a
            number_fixed = n
            number_unfixed = m
            assert (specie_a in atomic_symbols)
        else:
            atom_fixed = specie_b
            number_fixed = m
            number_unfixed = n
            assert (specie_b in atomic_symbols)

        ret = []
        for entry in self.entries.find({'structure.nspecies': 2, 'structure.formula': {'$regex': atom_fixed}}):
            comp = Structure.from_dict(entry['structure']).get_composition()
            if atom_fixed in comp.composition and comp.composition[atom_fixed] % number_fixed == 0:
                species = comp.species
                species.remove(atom_fixed)
                other_specie = species[0]
                # See if the fraction n/m is correct for A and B
                if number_unfixed * float(comp.composition[atom_fixed]) == number_fixed * float(
                        comp.composition[other_specie]):
                    ret.append(entry['_id'])
        return ret

    def find_composition(self, composition):
        """
        Search for structures with a pseudo-composition expressed as dictionary
        where symbols that are not atomic symbols such as A or X can be used to
        represent arbitrary atoms

        :return: (list) List of ids for all the structures that fulfill
                 the conditions
        """
        ret = []
        for entry in self.entries.find({'nspecies': len(composition)}):
            comp = Structure.from_dict(entry['structure']).get_composition()
            valid = True
            if sum(comp.composition.values()) % sum(composition.values()) != 0:
                valid = False
            vect1 = np.float64(np.sort(comp.composition.values()))
            vect2 = np.float64(np.sort(composition.values()))
            v12 = vect1 / vect2
            if not np.all(v12 == v12[0]):
                valid = False
            for ispecie in composition:
                if ispecie in atomic_symbols and ispecie not in comp.species:
                    valid = False

            if valid:
                print(comp.composition)
                ret.append(entry['_id'])
        return ret

    def get_structure(self, entry_id):
        """
        Return the structure in the entry with id 'entry_id'

        :rtype : Structure
        """
        entry = self.entries.find_one({'_id': entry_id})
        return Structure.from_dict(entry['structure'])

    def get_dicts(self, entry_id):
        """
        Return a tuple with the fields in an entry
        structure, properties and status

        :rtype : tuple
        """
        entry = self.entries.find_one({'_id': entry_id})
        structure_dict = entry['structure']
        if 'properties' in entry:
            properties = entry['properties']
        else:
            properties = None
        if 'status' in entry:
            status = entry['status']
        else:
            status = None

        return structure_dict, properties, status

    def get_entry(self, entry_id):
        return self.entries.find_one({'_id': entry_id})

    def is_locked(self, entry_id):
        """
        Return if a given entry is locked by someone
        evaluating the structure contained

        :rtype : bool
        """
        entry = self.entries.find_one({'_id': entry_id})

        if 'status' in entry and entry['status'] is not None and 'lock' in entry['status']:
            return True
        else:
            return False

    def lock(self, entry_id, name=None):
        if name is None:
            name = socket.gethostname()
        self.entries.update({'_id': entry_id}, {'$set': {'status.lock': name}})

    def unlock(self, entry_id, name=None):
        if name is None:
            self.entries.update({'_id': entry_id}, {'$unset': {'status.lock': 1}})
        else:
            self.entries.update({'_id': entry_id, 'status.lock': name}, {'$unset': {'status.lock': 1}})

    def set_minimal_schema(self):
        for entry_id in self.entries.find({'status': None}, {'status': 1}):
            print(entry_id)
            self.entries.update({'_id': entry_id['_id']}, {'$set': {'status': {}}})
        for entry_id in self.entries.find({'properties': None}):
            print(entry_id)
            self.entries.update({'_id': entry_id['_id']}, {'$set': {'properties': {}}})

    def create_static(self, field):

        for entry in self.entries.find({}):
            entry[field + '_static'] = entry[field]
            self.db.pychemia_entries.update({'_id': entry['_id']}, entry)

    def map_to_all(self, function, nparal=6):

        pool = Pool(processes=nparal)
        cursor = self.entries.find({}, no_cursor_timeout=True)
        print(cursor.count())
        entries = [entry['_id'] for entry in cursor]
        ret = pool.map(function, itertools.izip(itertools.repeat(self.db_settings), entries))
        cursor.close()
        return ret

    def replace_failed(self):
        for entry in self.entries.find({'status.relaxation': 'failed'}):
            st = self.get_structure(entry['_id'])
            comp = st.composition
            new_structure = Structure.random_cell(comp)
            self.entries.update({'_id': entry['_id']}, {'$unset': {'status.relaxation': 1,
                                                                   'status.target_forces': 1,
                                                                   'properties.energy': 1,
                                                                   'properties.forces': 1,
                                                                   'properties.stress': 1}})
            self.update(entry['_id'], structure=new_structure)

    def get_tags(self):
        entries = [x['status'].keys() for x in self.entries.find({}, {'status': 1})]
        ret = []
        for i in entries:
            for j in i:
                if j not in [u'target_forces', u'relaxation', u'lock'] and j not in ret:
                    ret.append(j)
        return ret


def get_database(db_settings):
    """
    Return a PyChemiaDB object either by recovering the database from its name on MongoDB or
    by creating a new one. The argument is a single python dictionary that should contain
    keys and values to create or get the database.

    """
    if 'host' not in db_settings:
        db_settings['host'] = 'localhost'
    if 'port' not in db_settings:
        db_settings['port'] = 27017
    if 'ssl' not in db_settings:
        db_settings['ssl'] = False
    if 'replicaset' not in db_settings:
        db_settings['replicaset'] = None
    if 'user' not in db_settings:
        pcdb = PyChemiaDB(name=db_settings['name'], host=db_settings['host'], port=db_settings['port'],
                          ssl=db_settings['ssl'], replicaset=db_settings['replicaset'])
    else:
        if 'admin_name' in db_settings and 'admin_passwd' in db_settings:
            pcdb = create_database(name=db_settings['name'], host=db_settings['host'], 
                                   port=db_settings['port'], ssl=db_settings['ssl'], 
                                   user_name=db_settings['user'], 
                                   user_passwd=db_settings['passwd'],
                                   admin_name=db_settings['admin_name'], 
                                   admin_passwd=db_settings['admin_passwd'],
                                   replicaset=db_settings['replicaset'])
        else:
            pcdb = PyChemiaDB(name=db_settings['name'], host=db_settings['host'], 
                              port=db_settings['port'], user=db_settings['user'], 
                              passwd=db_settings['passwd'], ssl=db_settings['ssl'],
                              replicaset=db_settings['replicaset'])
    return pcdb


def create_database(name, admin_name, admin_passwd, user_name, user_passwd, 
                    host='localhost', port=27017, ssl=False, replicaset=None):
    """
    Creates a new database for the database 'name'

    :param name: (str) The name of the database
    :param admin_name: (str) The administrator name
    :param admin_passwd: (str) Administrator password
    :param user_name: (str) Username for the database
    :param user_passwd: (str) Password for the user
    :param host: (str) Name of the host for the MongoDB server (default: 'localhost')
    :param port: (int) Port to connect to the MongoDB server (default: 27017)
    :param ssl: (bool) If True enable ssl encryption for communications to host (default: False)
    :param replicaset: (str, None) Identifier of a Replica Set

    """
    maxSevSelDelay = 2
    mc = pymongo.MongoClient(host=host, port=port, ssl=ssl, 
                             ssl_cert_reqs=pymongo.ssl_support.ssl.CERT_NONE,
                             replicaset=replicaset, serverSelectionTimeoutMS=maxSevSelDelay)
    try:
        mc.admin.authenticate(admin_name, admin_passwd)
    except OperationFailure:
        raise RuntimeError("Could not authenticate with the provided credentials")

    mc[name].add_user(user_name, user_passwd)
    return PyChemiaDB(name=name, user=user_name, passwd=user_passwd, host=host, port=port, ssl=ssl,
                      replicaset=replicaset)


def object_id(entry_id):
    if isinstance(entry_id, str):
        return ObjectId(entry_id)
    elif isinstance(entry_id, ObjectId):
        return entry_id


def has_connection(host='localhost'):
    if not HAS_PYMONGO:
        return False
    import pymongo
    try:
        maxSevSelDelay = 2
        client = pymongo.MongoClient(host, serverSelectionTimeoutMS=maxSevSelDelay)
        client.server_info()  # force connection on a request as the
        # connect=True parameter of MongoClient seems
        # to be useless here
        return True
    except pymongo.errors.ServerSelectionTimeoutError as err:
        # do whatever you need
        print(err)
        return False
