__author__ = 'Guillermo Avendano-Franco'

import os
import pychemia
import random
import numpy as np
from _population import Population
from pychemia import pcm_log


class PopulationNonColl(Population):

    def __init__(self, name, source_dir='.', mag_atoms=None):
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

    def __str__(self):
        ret = ' Population NonColl\n\n'
        ret += ' Name:               %s\n' % self.name
        ret += ' Tag:                %s\n' % self.tag
        ret += ' Formula:            %s\n' % self.structure.formula

        ret += ' Members:            %d\n' % len(self.members)
        ret += ' Actives:            %d\n' % len(self.actives)
        ret += ' Evaluated:          %d\n' % len(self.evaluated)
        return ret

    def add_random(self):
        """

        :return:
        """
        magmom = np.zeros((self.structure.natom, 3))
        mag_angles = []
        for i in range(self.structure.natom):
            if i in self.mag_atoms:
                angles = 2 * np.pi * np.random.random(3)
                mag = 2
                mag_angles.append([mag] + list(angles))
                rx = pychemia.utils.mathematics.rotation_x(angles[0])
                ry = pychemia.utils.mathematics.rotation_y(angles[1])
                rz = pychemia.utils.mathematics.rotation_z(angles[2])
                vec = mag * np.dot(rz, np.dot(ry, np.dot(rx, np.array([1, 0, 0]))))
                magmom[i] = vec
        data = {'mag_angles': mag_angles, 'magmom': list(magmom.flatten())}

        return self.new_entry(data)

    def cross(self, ids):
        pass

    def from_dict(self, population_dict):
        pass

    def new_entry(self, data, active=True):
        properties = {'mag_angles': data['mag_angles'], 'magmom': data['magmom']}
        status = {self.tag: active}
        entry_id = self.pcdb.insert(structure=self.structure, properties=properties, status=status)
        pcm_log.debug('Added new entry: %s with tag=%s: %s' % (str(entry_id), self.tag, str(active)))
        return entry_id

    def is_evaluated(self, entry_id):
        pass

    def check_duplicates(self, ids):
        pass

    def distance(self, entry_id, entry_jd):
        pass

    def move_random(self, entry_id, factor=0.2, in_place=False, kind='move'):
        pass

    def move(self, entry_id, entry_jd, factor=0.2, in_place=False):
        pass

    def recover(self):
        pass

    def value(self, entry_id):
        pass

    def str_entry(self, entry_id):
        entry = self.get_entry(entry_id)
        print entry['properties']['mag_angles']

    def get_duplicates(self, ids):
        return None
