import os
import numpy as np
from ._population import Population
from pychemia import pcm_log
from pychemia.utils.mathematics import spherical_to_cartesian, cartesian_to_spherical, rotate_towards_axis, \
    angle_between_vectors
from pychemia.code.vasp import read_incar, read_poscar, VaspJob, VaspOutput
from pychemia.crystal import KPoints


class NonCollinearMagMoms(Population):

    def evaluate_entry(self, entry_id):
        pass

    def __init__(self, name, source_dir='.', mag_atoms=None, magmom_magnitude=2.0, distance_tolerance=0.1, debug=False):
        Population.__init__(self, name, 'global')
        if not os.path.isfile(source_dir + os.sep + 'INCAR'):
            raise ValueError("INCAR not found")
        if not os.path.isfile(source_dir + os.sep + 'POSCAR'):
            raise ValueError("POSCAR not found")
        if not os.path.isfile(source_dir + os.sep + 'POTCAR'):
            raise ValueError("POTCAR not found")
        if not os.path.isfile(source_dir + os.sep + 'KPOINTS'):
            raise ValueError("KPOINTS not found")
        self.input = read_incar(source_dir + os.sep + 'INCAR')
        self.source_dir = source_dir

        if 'MAGMOM' not in self.input:
            raise ValueError('INCAR should define the MAGMOM variable')
        magmom = np.array(self.input.MAGMOM).reshape((-1, 3))

        self.structure = read_poscar(source_dir + os.sep + 'POSCAR')
        if mag_atoms is None:
            self.mag_atoms = list(np.where(np.apply_along_axis(np.linalg.norm, 1, magmom) > 0.0)[0])
            self.mag_atoms = [int(x) for x in self.mag_atoms]
        else:
            self.mag_atoms = mag_atoms
        self.magmom_magnitude = magmom_magnitude
        self.distance_tolerance = distance_tolerance
        self.debug = debug

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
    def from_dict(self, population_dict):
        return NonCollinearMagMoms(name=population_dict['name'],
                                 mag_atoms=population_dict['mag_atoms'],
                                 magmom_magnitude=population_dict['magmom_magnitude'],
                                 distance_tolerance=population_dict['distance_tolerance'])

    def debug_evaluation(self, magmom_sph):
        if self.debug:
            magmom_car=spherical_to_cartesian(magmom_sph)
            good_magmom = np.zeros((self.structure.natom, 3))
            for i in self.mag_atoms:
                good_magmom[i] = (-1**i) * 1.15470054 * np.ones(3)
            distance = np.sum(angle_between_vectors(magmom_car, good_magmom))
            distance /= len(self.mag_atoms)
            return distance-np.pi
        else:
            return None

    def new_entry(self, data, active=True):
        # Magnetic moments are stored in spherical coordinates
        data = np.array(data)
        properties = {'magmom': list(data.flatten()), 'energy': self.debug_evaluation(data)}
        status = {self.tag: active}
        entry={'structure': self.structure.to_dict, 'properties': properties, 'status': status}
        entry_id = self.insert_entry(entry)
        #pcm_log.debug('Added new entry: %s with tag=%s: %s' % (str(entry_id), self.tag, str(active)))
        return entry_id

    def is_evaluated(self, entry_id):
        entry = self.get_entry(entry_id, {'_id': 0, 'properties': 1})
        if 'energy' in entry['properties']:
            return True
        else:
            return False

    def check_duplicates(self, ids):
        selection = self.ids_sorted(ids)
        ret = {}
        for i in range(len(ids) - 1):
            for j in range(i+1, len(ids)):
                distance=self.distance(selection[i], selection[j])
                #print('The distance between candidates [%d, %d] is: %f' % (i,j,distance))
                if distance < self.distance_tolerance:
                    ret[selection[j]] = selection[i]
        return ret

    def distance(self, entry_id, entry_jd):
        entry = self.get_entry(entry_id, {'properties.magmom': 1})
        magmom_i = np.array(entry['properties']['magmom']).reshape((-1, 3))
        entry = self.get_entry(entry_jd, {'properties.magmom': 1})
        magmom_j = np.array(entry['properties']['magmom']).reshape((-1, 3))

        #print(entry_id)
        #print(magmom_i[self.mag_atoms])
        #print(entry_jd)
        #print(magmom_j[self.mag_atoms])

        magmom_ixyz = spherical_to_cartesian(magmom_i)
        magmom_jxyz = spherical_to_cartesian(magmom_j)
        distance = np.sum(angle_between_vectors(magmom_ixyz, magmom_jxyz))
        distance /= len(self.mag_atoms)
        return distance

    def move_random(self, entry_id, factor=0.2, in_place=False, kind='move'):

        entry = self.get_entry(entry_id, {'properties.magmom': 1})
        # Magnetic Momenta are stored in spherical coordinates
        magmom_i = spherical_to_cartesian(entry['properties']['magmom'])
        # Converted into cartesians
        magmom_xyz = spherical_to_cartesian(magmom_i)
        # Randomly disturbed using the factor
        magmom_xyz += factor * np.random.random((self.structure.natom, 3)) - factor / 2
        # Reconverting to spherical coordinates
        magmom_new = cartesian_to_spherical(magmom_xyz)
        # Resetting magnitudes
        magmom_new[:, 0] = self.magmom_magnitude

        properties = {'magmom': list(magmom_new.flatten()), 'energy': self.debug_evaluation(magmom_new)}

        if in_place:
            return self.update_properties(entry_id, new_properties=properties)
        else:
            return self.new_entry(magmom_new, active=False)

    def move(self, entry_id, entry_jd, factor=0.2, in_place=False):

#FAKE FACTOR
#        factor=1.0
        magmom_new_xyz = np.zeros((self.structure.natom, 3)) 
        entry = self.get_entry(entry_id, {'properties.magmom': 1})
        magmom_i = np.array(entry['properties']['magmom']).reshape((-1, 3))
        magmom_ixyz = spherical_to_cartesian(magmom_i)
        entry = self.get_entry(entry_jd, {'properties.magmom': 1})
        magmom_j = np.array(entry['properties']['magmom']).reshape((-1, 3))
        magmom_jxyz = spherical_to_cartesian(magmom_j)

#        print('Rotating %s towards %s with factor %f' % (entry_id, entry_jd, factor))
#        print('Origin:\n %s ' % magmom_i[self.mag_atoms])
#        print('Destin:\n %s ' % magmom_j[self.mag_atoms])
        for i in self.mag_atoms:
            #print('Atom %d' % i)
            if magmom_i[i][0] > 0 and magmom_j[i][0] > 0:         
                magmom_new_xyz[i] = rotate_towards_axis(magmom_ixyz[i], magmom_jxyz[i],
                                                                                   fraction=factor)
                #print('Final magmom Cartesian : %d     %s' % (i,magmom_new_xyz[i]))

        magmom_new = cartesian_to_spherical(magmom_new_xyz)
#        print('Final spherical:\n %s' % magmom_new[self.mag_atoms])
        magmom_new[:, 0] = self.magmom_magnitude

        properties = {'magmom': list(magmom_new.flatten()), 'energy': self.debug_evaluation(magmom_new)}

        if in_place:
            self.update_properties(entry_id, new_properties=properties)
            return entry_id
        else:
            return self.new_entry(magmom_new, active=False)

    def value(self, entry_id):
        entry = self.get_entry(entry_id, {'properties.energy': 1})
        if 'energy' in entry['properties']:
            return entry['properties']['energy']
        else:
            return None

    def str_entry(self, entry_id):
        entry = self.get_entry(entry_id, {'properties': 1})
        ret='['
        for i in range(len(self.mag_atoms)):
            ret+= ("[%8.5f %8.5f %8.5f] " % tuple(np.array(entry['properties']['magmom']).reshape((-1, 3))[i]))
        ret+='] '
        ret+='Energy= %f' % entry['properties']['energy']
        return ret

    def get_duplicates(self, ids):
        return None

    def add_random(self):
        """

        :return:
        """
        n = self.structure.natom
        a = self.magmom_magnitude * np.ones(n)
        b = 2 * np.pi * np.random.random(n) - np.pi
        c = np.pi * np.random.random(n)
        magmom = np.vstack((a, b, c)).T
        for i in range(self.structure.natom):
            if i not in self.mag_atoms:
                magmom[i, :] = 0.0

        return self.new_entry(magmom), None

    def recover(self):
        data = self.get_population_info()
        if data is not None:
            self.mag_atoms = data['mag_atoms']
            self.distance_tolerance = data['distance_tolerance']
            self.name = data['name']
            self.magmom_magnitude = data['magmom_magnitude']

    def cross(self, ids):
        entry_id = ids[0]
        entry_jd = ids[1]
        entry = self.get_entry(entry_id, {'properties.magmom': 1})
        magmom_i = np.array(entry['properties']['magmom']).reshape((-1, 3))
        entry = self.get_entry(entry_jd, {'properties.magmom': 1})
        magmom_j = np.array(entry['properties']['magmom']).reshape((-1, 3))
        magmom_inew = np.zeros((self.structure.natom, 3))
        magmom_jnew = np.zeros((self.structure.natom, 3))
        for i in range(self.structure.natom):
            rnd = np.random.random()
            if rnd < 0.5:
                magmom_inew[i] = magmom_j[i]
                magmom_jnew[i] = magmom_i[i]
            else:
                magmom_inew[i] = magmom_i[i]
                magmom_jnew[i] = magmom_j[i]

        entry_id = self.new_entry(magmom_inew, active=True)
        entry_jd = self.new_entry(magmom_jnew, active=True)

        return entry_id, entry_jd

    def prepare_folder(self, entry_id, workdir, binary='vasp', source_dir='.'):

        for i in ['KPOINTS', 'POSCAR', 'POTCAR']:
            os.symlink(self.source_dir+os.sep+i , workdir+os.sep+i)

        input = read_incar(self.source_dir + os.sep + 'INCAR')
        magmom_sph = self.get_entry(entry_id, {'properties.magmom': 1})['properties']['magmom']
        magmom_car = spherical_to_cartesian(magmom_sph)
        input['MAGMOM'] = [float(x) for x in magmom_car.flatten()]
        input['M_CONSTR'] = [float(x) for x in magmom_car.flatten()]
        input['IBRION'] = -1
        input['LWAVE'] = True
        input['EDIFF'] = 1E-5
        input['LAMBDA'] = 10
        input['NSW'] = 0
        input['I_CONSTRAINED_M'] = 1
        input.write(workdir+ os.sep + 'INCAR')

    def collect_data(self, entry_id, workdir):
        if os.path.isfile(workdir + '/OUTCAR'):
            vo = VaspOutput(workdir + '/OUTCAR')
            if 'energy' in vo.final_data:
                if 'free_energy' in vo.final_data['energy']:
                    energy = vo.final_data['energy']['free_energy']
                    print('Uploading energy data for %s' % entry_id)
                    self.set_in_properties(entry_id, 'energy', energy)
                    return True
                else:
                    return False
            else:
                return False
        else:
            return False
