__author__ = 'Guillermo Avendano Franco'

from pymongo import MongoClient
from bson.objectid import ObjectId
import numpy as np
import socket

from pychemia.utils.periodic import atomic_symbols
from pychemia import Structure


class PyChemiaDB():

    def __init__(self, name='pychemiadb', host='localhost', port=27017, user=None, password=None):
        """
        Creates a MongoDB client to 'host' with 'port' and connect it to the database 'name'.
        Authentication can be used with 'user' and 'password'

        :param name: (str) The name of the database
        :param host: (str) The host as name or IP
        :param port: (int) The number of port to connect with the server (Default is 27017)
        :param user: (str) The user with read or write permissions to the database
        :param password: (str/int) Password to authenticate the user into the server

        :return:
        """
        self.name = name
        uri = 'mongodb://'
        if user is not None:
            uri += user
            if password is not None:
                uri += ':'+str(password)
            uri += '@'
        uri += host + ':' + str(port)
        if user is not None:
            uri += '/' + name
        self._client = MongoClient(uri)
        self.db = self._client[name]
        self.entries = self.db.pychemia_entries

    def insert(self, structure, properties=None, status=None):
        """
        Insert a pychemia structure instance and properties
        into the database
        :param structure: (pychemia.Structure) An instance of Pychemia's Structure
        :param properties: (dict) Dictionary of properties
        :param status: (dict) Dictionary of status
        :return:
        """
        entry = {'structure': structure.to_dict(), 'properties': properties, 'status': status}
        entry_id = self.entries.insert(entry)
        return entry_id

    def clean(self):
        self._client.drop_database(self.name)
        self.db = self._client[self.name]

    def update(self, entry_id, structure=None, properties=None, status=None):

        entry_id = object_id(entry_id)

        assert (self.entries.find_one({'_id': entry_id}) is not None)
        entry = self.entries.find_one({'_id': entry_id})
        if structure is not None:
            if isinstance(structure, Structure):
                entry['structure'] = structure.to_dict()
            elif isinstance(structure, dict):
                entry['structure'] = structure
            else:
                print 'ERROR: Could not process the structure'
                print type(structure)
        if properties is not None:
            entry['properties'] = properties
        if status is not None:
            entry['status'] = status

        self.entries.update({'_id': entry_id}, entry)

    def find_AnBm(self, specie_a=None, specie_b=None, n=1, m=1):
        """
        Search for structures with a composition expressed as AnBm
        where one and only one between A or B is fixed and the numbers
        amounts n and m are both fixed

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
        for entry in self.entries.find({'nspecies': 2, 'formula': {'$regex': atom_fixed}}):
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
                print comp.composition
                ret.append(entry['_id'])
        return ret

    def get_structure(self, entry_id):
        entry_id = object_id(entry_id)
        entry = self.entries.find_one({'_id': entry_id})
        return Structure.from_dict(entry['structure'])

    def get_dicts(self, entry_id):
        entry_id = object_id(entry_id)
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

    def is_locked(self, entry_id):

        entry_id = object_id(entry_id)
        entry = self.entries.find_one({'_id': entry_id})

        if 'status' in entry and 'lock' in entry['status']:
            return True
        else:
            return False

    def lock(self, entry_id):
        entry_id = object_id(entry_id)
        structure, properties, status = self.get_dicts(entry_id)
        status['lock'] = socket.gethostname()
        self.update(entry_id, status=status)

    def unlock(self, entry_id, name=None):
        entry_id = object_id(entry_id)
        structure, properties, status = self.get_dicts(entry_id)
        lockedby = None
        if 'lock' in status:
            if name is None:
                lockedby = status.pop('lock')
            elif status['lock'] == name:
                lockedby = status.pop('lock')
        self.update(entry_id, status=status)
        return lockedby


def get_database(db_settings):

    if 'host' not in db_settings:
        db_settings['host'] = 'localhost'
    if 'port' not in db_settings:
        db_settings['port'] = 27017

    if 'user' not in db_settings:
        db = PyChemiaDB(name=db_settings['name'], host=db_settings['host'], port=db_settings['port'])
    else:
        db = PyChemiaDB(name=db_settings['name'], host=db_settings['host'], port=db_settings['port'],
                        user=db_settings['user'], password=db_settings['password'])
    return db


def object_id(entry_id):
    if isinstance(entry_id, basestring):
        return ObjectId(entry_id)
    elif isinstance(entry_id, ObjectId):
        return entry_id