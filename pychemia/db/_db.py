__author__ = 'Guillermo Avendano Franco'

from pymongo import MongoClient
from bson.objectid import ObjectId
import numpy as np

from pychemia.utils.periodic import atomic_symbols
from pychemia import Structure


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
        entry_dict = structure.to_dict()
        if properties is not None:
            entry_dict['properties'] = properties
        else:
            entry_dict['properties'] = {}
        entry_id = self.entries.insert(entry_dict)
        return entry_id

    def get_iterator(self):
        cursor = self.db.structures.find()
        return Iterator(self.db, cursor)

    def clean(self):
        self._client.drop_database(self.name)
        self.db = self._client[self.name]

    @property
    def is_master(self):
        return self._is_master

    def update(self, entry_id, new_entry):

        if isinstance(entry_id, basestring):
            entry_id = ObjectId(entry_id)
        elif isinstance(entry_id, ObjectId):
            entry_id = entry_id

        assert (self.entries.find_one({'_id': entry_id}) is not None)
        self.entries.update({'_id': entry_id}, new_entry)

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
            comp = Structure.from_dict(entry).get_composition()
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
            comp = Structure.from_dict(entry).get_composition()
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
