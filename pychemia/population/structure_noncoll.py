import os

import numpy as np

import pychemia
from ._population import Population
from pychemia import pcm_log
from pychemia.utils.computing import deep_unicode

class PopulationNonColl(Population):
    def __init__(self, name, source_dir='.', mag_atoms=None, magmom_magnitude=2.0, distance_tolerance=0.1):
        Population.__init__(self, name, 'global')
        if not os.path.isfile(source_dir + os.sep + 'INCAR'):
            raise ValueError("INCAR not found")
        if not os.path.isfile(source_dir + os.sep + 'POSCAR'):
            raise ValueError("POSCAR not found")
        self.input = pychemia.code.vasp.read_incar(source_dir + os.sep + 'INCAR')
        magmom = np.array(self.input.get_value('MAGMOM')).reshape((-1, 3))

        self.structure = pychemia.code.vasp.read_poscar(source_dir + os.sep + 'POSCAR')
        if mag_atoms is None:
            self.mag_atoms = list(np.where(np.apply_along_axis(np.linalg.norm, 1, magmom) > 0.0)[0])
        else:
            self.mag_atoms = mag_atoms
        self.magmom_magnitude = magmom_magnitude
        self.distance_tolerance = distance_tolerance


    def __str__(self):
        ret = ' Population NonColl\n\n'
        ret += ' Name:               %s\n' % self.name
        ret += ' Tag:                %s\n' % self.tag
        ret += ' Formula:            %s\n' % self.structure.formula

        ret += ' Members:            %d\n' % len(self.members)
        ret += ' Actives:            %d\n' % len(self.actives)
        ret += ' Evaluated:          %d\n' % len(self.evaluated)
        return ret

    @property
    def to_dict(self):
        return {'name': self.name,
                'tag': self.tag,
                'mag_atoms': self.mag_atoms,
                'magmom_magnitude': self.magmom_magnitude,
                'distance_tolerance': self.distance_tolerance}

    @staticmethod
    def from_dict(population_dict):
        return PopulationNonColl(name=population_dict['name'],
                                 mag_atoms=population_dict['mag_atoms'],
                                 magmom_magnitude=population_dict['magmom_magnitude'],
                                 distance_tolerance=population_dict['distance_tolerance'])

    def new_entry(self, data, active=True):
        data = np.array(data)
        # Magnetic moments are stored in spherical coordinates
        properties = {'magmom': list(data.flatten())}
        status = {self.tag: active}
        entry_id = self.pcdb.insert(structure=self.structure, properties=properties, status=status)
        pcm_log.debug('Added new entry: %s with tag=%s: %s' % (str(entry_id), self.tag, str(active)))
        return entry_id

    def is_evaluated(self, entry_id):
        structure, properties, status = self.pcdb.get_dicts(entry_id)
        if 'energy' in properties:
            return True
        else:
            return False

    def check_duplicates(self, ids):
        selection = self.ids_sorted(ids)
        ret = {}
        for i in range(len(ids) - 1):
            for j in range(i, len(ids)):
                if self.distance(selection[i], selection[j]) < self.distance_tolerance:
                    ret[selection[j]] = selection[i]
        return ret

    def distance(self, entry_id, entry_jd):
        entry = self.pcdb.entries.find_one({'_id': entry_id}, {'properties.magmom': 1})
        magmom_i = pychemia.utils.mathematics.spherical_to_cartesian(entry['properties']['magmom'])
        entry = self.pcdb.entries.find_one({'_id': entry_jd}, {'properties.magmom': 1})
        magmom_j = pychemia.utils.mathematics.spherical_to_cartesian(entry['properties']['magmom'])
        magmom_ixyz = pychemia.utils.mathematics.spherical_to_cartesian(magmom_i)
        magmom_jxyz = pychemia.utils.mathematics.spherical_to_cartesian(magmom_j)
        distance = np.sum(pychemia.utils.mathematics.angle_between_vectors(magmom_ixyz, magmom_jxyz))
        distance /= len(self.mag_atoms)
        return distance

    def move_random(self, entry_id, factor=0.2, in_place=False, kind='move'):

        entry = self.pcdb.entries.find_one({'_id': entry_id}, {'properties.magmom': 1})
        # Magnetic Momenta are stored in spherical coordinates
        magmom_i = pychemia.utils.mathematics.spherical_to_cartesian(entry['properties']['magmom'])
        # Converted into cartesians
        magmom_xyz = pychemia.utils.mathematics.spherical_to_cartesian(magmom_i)
        # Randomly disturbed using the factor
        magmom_xyz += factor * np.random.rand((self.structure.natom, 3)) - factor / 2
        # Reconverting to spherical coordinates
        magmom_new = pychemia.utils.mathematics.cartesian_to_spherical(magmom_xyz)
        # Resetting magnitudes
        magmom_new[:, 0] = self.magmom_magnitude
        properties = {'magmom': magmom_new}

        if in_place:
            return self.pcdb.update(entry_id, properties=properties)
        else:
            return self.new_entry(magmom_new, active=False)

    def move(self, entry_id, entry_jd, factor=0.2, in_place=False):
        magmom_new_xyz = np.zeros((self.structure.natom, 3))
        entry = self.pcdb.entries.find_one({'_id': entry_id}, {'properties.magmom': 1})
        magmom_i = np.array(entry['properties']['magmom']).reshape((-1, 3))
        magmom_ixyz = pychemia.utils.mathematics.spherical_to_cartesian(magmom_i)
        entry = self.pcdb.entries.find_one({'_id': entry_jd}, {'properties.magmom': 1})
        magmom_j = np.array(entry['properties']['magmom']).reshape((-1, 3))
        magmom_jxyz = pychemia.utils.mathematics.spherical_to_cartesian(magmom_j)

        for i in range(self.structure.natom):
            if magmom_ixyz[i][0] > 0 and magmom_jxyz[i][0] > 0:
                magmom_new_xyz[i] = pychemia.utils.mathematics.rotate_towards_axis(magmom_ixyz[i], magmom_jxyz[i],
                                                                                   fraction=factor)

        magmom_new = pychemia.utils.mathematics.cartesian_to_spherical(magmom_new_xyz)
        magmom_new[:, 0] = self.magmom_magnitude
        properties = {'magmom': magmom_new}

        if in_place:
            return self.pcdb.update(entry_id, properties=properties)
        else:
            return self.new_entry(magmom_new, active=False)

    def value(self, entry_id):
        entry = self.pcdb.entries.find_one({'_id': entry_id}, {'properties.energy': 1})
        if 'energy' in entry['properties']:
            return entry['properties']['energy']
        else:
            return None

    def str_entry(self, entry_id):
        entry = self.pcdb.entries.find_one({'_id': entry_id}, {'properties.magmom': 1})
        print(np.array(entry['properties']['magmom']).reshape((-1, 3)))

    def get_duplicates(self, ids):
        return None

    def add_random(self):
        """

        :return:
        """
        n = self.structure.natom
        a = self.magmom_magnitude * np.ones(n)
        b = 2 * np.pi * np.random.rand(n)
        c = np.pi * np.random.rand(20)
        magmom = np.vstack((a, b, c)).T
        for i in range(self.structure.natom):
            if i not in self.mag_atoms:
                magmom[i, :] = 0.0

        return self.new_entry(magmom)

    def recover(self):
        data = self.pcdb.db.population_info.find_one({'tag': self.tag})
        if data is not None:
            self.mag_atoms = data['mag_atoms']
            self.distance_tolerance = data['distance_tolerance']
            self.name = data['name']
            self.magmom_magnitude = data['magmom_magnitude']

    def cross(self, ids):
        entry_id = ids[0]
        entry_jd = ids[1]
        entry = self.pcdb.entries.find_one({'_id': entry_id}, {'properties.magmom': 1})
        magmom_i = np.array(entry['properties']['magmom']).reshape((-1, 3))
        entry = self.pcdb.entries.find_one({'_id': entry_jd}, {'properties.magmom': 1})
        magmom_j = np.array(entry['properties']['magmom']).reshape((-1, 3))
        magmom_inew = np.zeros((self.structure.natom, 3))
        magmom_jnew = np.zeros((self.structure.natom, 3))
        for i in range(self.structure.natom):
            rnd = np.random.rand()
            if rnd < 0.5:
                magmom_inew[i] = magmom_j[i]
                magmom_jnew[i] = magmom_i[i]
            else:
                magmom_inew[i] = magmom_i[i]
                magmom_jnew[i] = magmom_j[i]

        entry_id = self.new_entry(magmom_inew, active=True)
        entry_jd = self.new_entry(magmom_jnew, active=True)

        return entry_id, entry_jd
