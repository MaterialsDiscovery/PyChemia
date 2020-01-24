import os
import numpy as np
from ._population import Population
from pychemia.utils.mathematics import spherical_to_cartesian, cartesian_to_spherical, rotate_towards_axis, \
    angle_between_vectors
from pychemia.code.vasp import read_incar, read_poscar, VaspOutput


class NonCollinearMagMoms(Population):

    def __init__(self, name, source_dir='.', mag_atoms=None, magmom_magnitude=2.0, distance_tolerance=0.1,
                 incar_extra=None, debug=False):
        """
        This class provides a population of Magnetic Moment vectors for the same structure and was created to be used
        on VASP.
        The magnetic moment is set in cartesian coordinates with 3 numbers for each atom in the unit cell.
        This population provides methods to manipulate the magnetic moments between different candidates in order to
        optimize the magnetic orientations using the global-search methods implemented on PyChemia.

        :param name: The name of the database to be created or directly the PyChemiaDB database object.
                     The name is used when the database can be created without username, password and no encryption.
                     Otherwise the database must be created first and its object be 'name' argument.
        :param source_dir: Directory that contains the basic 4 files for VASP: 'POSCAR', 'POTCAR', 'KPOINTS' and 'INCAR'
                            Except for 'INCAR' the files are linked symbolically on each directory that will run VASP
                            The input variables on 'INCAR' changing only MAGMOM and I_CONSTRAINED_M.
                            The 'INCAR' file could contain generic variables and some other variables could be directly
                            specified using the dictionary 'incar_extra'.
        :param mag_atoms: List of atoms for which the Magnetic Moments are changed. If the variable is None, the list
                            is inferred from the original INCAR file. The numbering of atoms start with 0.
        :param magmom_magnitude: Fix value for the magnitude of the magnetic moment imposed for all the atoms in
                                'mag_atoms' list
        :param distance_tolerance: Maximal distance in magnetic moments to consider two candidates as equivalent.
        :param incar_extra: Extra variables for INCAR file that are added or replaced from the INCAR read with
                            'source_dir'
        :param debug: If True produce a verbose output during the different calls to the methods.
        """
        Population.__init__(self, name, 'global', distance_tolerance=distance_tolerance)
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
        if incar_extra is None:
            self.incar_extra = {'IBRION': -1,
                                'LWAVE': True,
                                'LAMBDA': 10,
                                'NSW': 0,
                                'I_CONSTRAINED_M': 1}

        else:
            self.incar_extra = incar_extra
        self.debug = debug

    def __str__(self):
        ret = ' Population NonColl\n\n'
        ret += ' Name:               %s\n' % self.name
        ret += ' Tag:                %s\n' % self.tag
        ret += ' Formula:            %s\n' % self.structure.formula
        ret += ' source_dir:         %s\n' % self.source_dir
        ret += ' magmom_magnitude:   %s\n' % self.magmom_magnitude
        ret += ' distance_tolerance: %s\n' % self.distance_tolerance
        ret += ' incar_extra:        %s\n' % self.incar_extra

        ret += ' Members:            %d\n' % len(self.members)
        ret += ' Actives:            %d\n' % len(self.actives)
        ret += ' Evaluated:          %d\n' % len(self.evaluated)
        return ret

    @property
    def to_dict(self):
        return {'name': self.name,
                'tag': self.tag,
                'mag_atoms': self.mag_atoms,
                'source_dir': self.source_dir,
                'incar_extra': self.incar_extra,
                'magmom_magnitude': self.magmom_magnitude,
                'distance_tolerance': self.distance_tolerance}

    @staticmethod
    def from_dict(self, population_dict):
        return NonCollinearMagMoms(name=population_dict['name'],
                                   source_dir=population_dict['source_dir'],
                                   mag_atoms=population_dict['mag_atoms'],
                                   incar_extra=population_dict['incar_extra'],
                                   magmom_magnitude=population_dict['magmom_magnitude'],
                                   distance_tolerance=population_dict['distance_tolerance'])

    def debug_evaluation(self, magmom_sph):
        """
        For debugging ONLY:
        Fake evaluation of total energy using a magnetic configuration as the minimal
        energy, the numerical distance with other candidates defines their energy.

        :param magmom_sph:
        :return:
        """
        if self.debug:
            magmom_car = spherical_to_cartesian(magmom_sph)
            good_magmom = np.zeros((self.structure.natom, 3))
            for i in self.mag_atoms:
                good_magmom[i] = (-1**i) * 1.15470054 * np.ones(3)
            distance = np.sum(angle_between_vectors(magmom_car, good_magmom))
            distance /= len(self.mag_atoms)
            return distance-np.pi
        else:
            return None

    def evaluate_entry(self, entry_id):
        pass

    def new_entry(self, cartesian_magmoms, active=True):
        """
        Creates a new entry on the database

        :param cartesian_magmoms: Magnetic moments stored in spherical coordinates
        :param active:
        :return:
        """
        cartesian_magmoms = np.array(cartesian_magmoms)
        properties = {'magmom': list(cartesian_magmoms.flatten()), 'energy': self.debug_evaluation(cartesian_magmoms)}
        status = {self.tag: active}
        entry = {'structure': self.structure.to_dict, 'properties': properties, 'status': status}
        entry_id = self.insert_entry(entry)
        return entry_id

    def is_evaluated(self, entry_id):
        """
        One candidate is considered evaluated if it contains any finite value of energy on the properties.energy field

        :param entry_id:
        :return:
        """
        entry = self.get_entry(entry_id, {'_id': 0, 'properties': 1})
        if 'energy' in entry['properties'] and entry['properties']['energy'] is not None:
            return True
        else:
            return False

    def distance(self, entry_id, entry_jd):
        """
        Compute the distance between 2 candidates as the average angular movement
        between the magnetic moments for all the atoms considered magnetic ('mag_atoms')

        :param entry_id:
        :param entry_jd:
        :return:
        """
        entry = self.get_entry(entry_id, {'properties.magmom': 1})
        magmom_i = np.array(entry['properties']['magmom']).reshape((-1, 3))
        entry = self.get_entry(entry_jd, {'properties.magmom': 1})
        magmom_j = np.array(entry['properties']['magmom']).reshape((-1, 3))

        magmom_ixyz = spherical_to_cartesian(magmom_i)
        magmom_jxyz = spherical_to_cartesian(magmom_j)
        distance = np.sum(angle_between_vectors(magmom_ixyz, magmom_jxyz))
        distance /= len(self.mag_atoms)
        return distance

    def move_random(self, entry_id, factor=0.2, in_place=False, kind='move'):
        """
        Change the magnetic orientation randomly for all the atoms considered magnetic 'mag_atoms'
        The 'factor' argument scales the intensite of the movement.

        :param entry_id:
        :param factor:
        :param in_place:
        :param kind:
        :return:
        """
        entry = self.get_entry(entry_id, {'properties.magmom': 1})
        # Magnetic Momenta are stored in spherical coordinates
        magmom_xyz = spherical_to_cartesian(entry['properties']['magmom'])
        # Randomly disturbed using the factor
        magmom_xyz += factor * np.random.random((self.structure.natom, 3)) - factor / 2
        # Reconverting to spherical coordinates
        magmom_new = cartesian_to_spherical(magmom_xyz)
        # Resetting magnitudes
        for i in range(len(magmom_xyz)):
            if magmom_xyz[i][0] > 0.0:
                magmom_new[i, 0] = self.magmom_magnitude
            else:
                magmom_new[i, 0] = 0.0

        properties = {'magmom': list(magmom_new.flatten()), 'energy': self.debug_evaluation(magmom_new)}

        if in_place:
            return self.update_properties(entry_id, new_properties=properties)
        else:
            return self.new_entry(magmom_new, active=False)

    def move(self, entry_id, entry_jd, factor=0.2, in_place=False):
        """
        Move the magnetic moments from one candidate in the direction of another.

        :param entry_id: Source candidate, this one will be moved
        :param entry_jd: Destination candidate, this one is not moved.
        :param factor: (float) a value in range [0,1] 0 being the source, 1 the destination
        :param in_place: (bool) if true the new magnetic moments replace those in entry_id
        :return:
        """
        magmom_new_xyz = np.zeros((self.structure.natom, 3))
        entry = self.get_entry(entry_id, {'properties.magmom': 1})
        magmom_i = np.array(entry['properties']['magmom']).reshape((-1, 3))
        magmom_ixyz = spherical_to_cartesian(magmom_i)
        entry = self.get_entry(entry_jd, {'properties.magmom': 1})
        magmom_j = np.array(entry['properties']['magmom']).reshape((-1, 3))
        magmom_jxyz = spherical_to_cartesian(magmom_j)

        for i in self.mag_atoms:
            # print('Atom %d' % i)
            if magmom_i[i][0] > 0 and magmom_j[i][0] > 0:
                magmom_new_xyz[i] = rotate_towards_axis(magmom_ixyz[i], magmom_jxyz[i], fraction=factor)
                # print('Final magmom Cartesian : %d     %s' % (i,magmom_new_xyz[i]))

        magmom_new = cartesian_to_spherical(magmom_new_xyz)

        magmom_new[:, 0] = self.magmom_magnitude
        for i in range(self.structure.natom):
            if i not in self.mag_atoms:
                magmom_new[i, :] = 0.0

        properties = {'magmom': list(magmom_new.flatten()), 'energy': self.debug_evaluation(magmom_new)}

        if in_place:
            self.update_properties(entry_id, new_properties=properties)
            return entry_id
        else:
            return self.new_entry(magmom_new, active=False)

    def value(self, entry_id):
        """
        Return the energy value associated to the candidate with identifier 'entry_id'

        :param entry_id:
        :return:
        """
        entry = self.get_entry(entry_id, {'properties.energy': 1})
        if 'energy' in entry['properties']:
            return entry['properties']['energy']
        else:
            return None

    def str_entry(self, entry_id):
        entry = self.get_entry(entry_id, {'properties': 1})
        ret = '['
        for i in range(len(self.mag_atoms)):
            ret += ("[%8.5f %8.5f %8.5f] " % tuple(np.array(entry['properties']['magmom']).reshape((-1, 3))[i]))
        ret += '] '
        ret += 'Energy= %f' % entry['properties']['energy']
        return ret

    def add_random(self):
        """
        Creates a new candidate with random orientations for the magnetic moments of all the atoms in 'mag_atoms'

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
            self.source_dir = data['source_dir']
            self.incar_extra = data['incar_extra']

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

    def prepare_folder(self, entry_id, workdir, executable='vasp', source_dir='.'):

        if not os.path.isdir(workdir):
            os.mkdir(workdir)

        for i in ['KPOINTS', 'POSCAR', 'POTCAR']:
            if os.path.exists(workdir+os.sep+i):
                os.remove(workdir+os.sep+i)
            os.symlink(os.path.abspath(self.source_dir+os.sep+i), workdir+os.sep+i)

        incar = read_incar(self.source_dir + os.sep + 'INCAR')
        magmom_sph = self.get_entry(entry_id, {'properties.magmom': 1})['properties']['magmom']
        magmom_car = spherical_to_cartesian(magmom_sph)
        incar['MAGMOM'] = [float(x) for x in magmom_car.flatten()]
        incar['M_CONSTR'] = [float(x) for x in magmom_car.flatten()]
        for i in self.incar_extra:
            incar[i] = self.incar_extra[i]
        incar.write(workdir + os.sep + 'INCAR')

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
