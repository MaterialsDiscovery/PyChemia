from __future__ import print_function
import os
import re
import itertools
import numpy as np
from ._population import Population
from pychemia import pcm_log
from pychemia.code.abinit import InputVariables
from pychemia.code.abinit import AbinitOutput
from pychemia.utils.mathematics import gram_smith_qr, gea_all_angles, gea_orthogonal_from_angles


class OrbitalDFTU(Population):
    def __init__(self, name, abinit_input='abinit.in', num_electrons_dftu=None, num_indep_matrices=None,
                 connections=None):

        """
        This population is created with the purpose of global optimization of correlated orbitals 'd'
        and 'f'.


        The population is basically a collection of ABINIT inputs, the candidates have the same structure and
        uses the same input variables with exception of 'dmatpawu', the purpose of the
        population is to use global-population searchers to find the correlation matrices
        'dmatpawu' that minimizes the energy.

        The variable 'dmatpawu' is a list of numbers that can be arranged into N matrices.
        The matrices are 5x5 for 'd' orbitals and 7x7 for 'f' orbitals.

        :param name: The name of the 'PyChemiaDB' database created to stored the different
                        set of variables and the resulting output from abinit.
                     When using databases protected with username and password, the database
                        should be created independently and the database object must be use
                        as the 'name' argument

        :param abinit_input: The abinit input file, all the variables will be preserve for all new candidates, except
                        for 'dmatpawu' the only variable that changes.

        :param num_electrons_dftu: Example [5, 1, 5, 1]

        """
        # Call the parent class initializer to link the PychemiaDB that will be used
        Population.__init__(self, name, 'global')
        # Checking for existence of 'abinit.in'
        if not os.path.isfile(abinit_input):
            print("Abinit input not found")
            raise ValueError
        # Reading the input variables and getting the structure
        self.input = InputVariables(abinit_input)
        self.structure = self.input.get_structure()

        print('Orbital population:')
        print('Species [znucl]: %s' % self.input['znucl'])

        self.natpawu = 0
        print('Orbitals corrected:')
        for i in range(self.input['ntypat']):
            if self.input['lpawu'][i] == -1:
                print("%3s : False" % self.input['znucl'][i])
            else:
                print("%3s : True (l=%d)" % (self.input['znucl'][i], self.input['lpawu'][i]))
                self.natpawu += sum([1 for x in self.input['typat'] if x == i + 1])
        print('Number of atoms where DFT+U is applied: %d' % self.natpawu)

        # Computing the orbital that will be corrected
        # 2 (d-orbitals) or 3 (f-orbitals)
        self.maxlpawu = max(self.input['lpawu'])
        if self.maxlpawu == 2:
            print("Correlation of 'd' orbitals")
        elif self.maxlpawu == 3:
            print("Correlation of 'f' orbitals")

        # nsppol is the number of independent spin polarisations. Can take the values 1 or 2
        if self.input.has_variable('nsppol'):
            self.nsppol = self.input.get_value('nsppol')
        else:
            # Default from ABINIT
            self.nsppol = 1

        # nspinor it the umber of spinorial components of the wavefunctions
        if self.input.has_variable('nspinor'):
            self.nspinor = self.input.get_value('nspinor')
        else:
            self.nspinor = 1

        # nspden is the number of spin-density components
        if self.input.has_variable('nspden'):
            self.nspden = self.input.get_value('nspden')
        else:
            self.nspden = self.nsppol

        if self.nsppol == 1 and self.nspinor == 1 and self.nspden == 1:
            # Non-magnetic system (nsppol=1, nspinor=1, nspden=1):
            # One (2lpawu+1)x(2lpawu+1) dmatpawu matrix is given for each atom on which +U is applied.
            # It contains the "spin-up" occupations.
            self.nmatrices = self.natpawu
        elif self.nsppol == 2 and self.nspinor == 1 and self.nspden == 2:
            # Ferromagnetic spin-polarized (collinear) system (nsppol=2, nspinor=1, nspden=2):
            # Two (2lpawu+1)x(2lpawu+1) dmatpawu matrices are given for each atom on which +U is applied.
            # They contain the "spin-up" and "spin-down" occupations.
            self.nmatrices = 2 * self.natpawu
        elif self.nsppol == 1 and self.nspinor == 1 and self.nspden == 2:
            # Anti-ferromagnetic spin-polarized(collinear) system(nsppol=1, nspinor=1, nspden=2):
            # One(2lpawu + 1)x(2lpawu + 1) dmatpawu matrix is given for each atom on which +U is applied.
            # It contains the "spin-up" occupations.
            self.nmatrices = self.natpawu
        elif self.nsppol == 1 and self.nspinor == 2 and self.nspden == 4:
            # Non-collinear magnetic system (nsppol=1, nspinor=2, nspden=4):
            # Two (2lpawu+1)x(2lpawu+1) dmatpawu matrices are given for each atom on which +U is applied.
            # They contains the "spin-up" and "spin-down" occupations (defined as n_up=(n+|m|)/2 and n_dn=(n-|m|)/2),
            #    where m is the integrated magnetization vector).
            self.nmatrices = 2 * self.natpawu
        elif self.nsppol == 1 and self.nspinor == 2 and self.nspden == 1:
            # Non-collinear magnetic system with zero magnetization (nsppol=1, nspinor=2, nspden=1):
            # Two (2lpawu+1)x(2lpawu+1) dmatpawu matrices are given for each atom on which +U is applied.
            # They contain the "spin-up" and "spin-down" occupations;
            self.nmatrices = 2 * self.natpawu

        print('Variables controling the total number of matrices')
        print('nsppol : %d' % self.nsppol)
        print('nspinor: %d' % self.nspinor)
        print('nspden : %d' % self.nspden)
        print('Total number of matrices expected on dmatpawu: %d' % self.nmatrices)

        if num_electrons_dftu is None:
            abiinput = InputVariables(abinit_input)
            dmatpawu = np.array(abiinput['dmatpawu']).reshape(-1, self.ndim, self.ndim)
            params = dmatpawu2params(dmatpawu, 5)
            self.num_electrons_dftu = np.apply_along_axis(sum,1,params['occupations'])
        else:
            self.num_electrons_dftu = np.array(num_electrons_dftu)
        print('Number of electrons for each correlation matrix: %s' % self.num_electrons_dftu)

        if num_indep_matrices is not None:
            self.num_indep_matrices = num_indep_matrices
        else:
            self.num_indep_matrices = self.nmatrices
        print('Number of independent matrices: %d' % self.num_indep_matrices)

        if connections is not None:
            self.connections = list(connections)
            if len(self.connections) != self.nmatrices:
                raise ValueError('Number of connections between matrices is not consistent with the number of matrices '
                                 'defined on dmatpawu')
            print('Connections: %s' % self.connections)
        else:
            self.connections = list(range(self.nmatrices))

    def __str__(self):
        ret = ' Population LDA+U\n\n'
        ret += ' Name:               %s\n' % self.name
        ret += ' Tag:                %s\n' % self.tag
        ret += ' Formula:            %s\n' % self.structure.formula
        ret += ' natpawu:            %d\n' % self.natpawu
        ret += ' nmatrices:          %d\n' % self.nmatrices
        ret += ' maxlpawu:           %d\n' % self.maxlpawu
        ret += ' num_indep_matrices: %s\n' % self.num_indep_matrices
        ret += ' connections:        %s\n' % self.connections

        ret += ' Members:            %d\n' % len(self.members)
        ret += ' Actives:            %d\n' % len(self.actives)
        ret += ' Evaluated:          %d\n' % len(self.evaluated)
        return ret

    @property
    def ndim(self):
        """
        Dimension of the matrices defined on dmatpawu, for 'd' orbitals is 5 for 'f' orbitals is 7

        :return:
        """
        return 2 * self.maxlpawu + 1

    def add_random(self):
        """
        Creates a new set of variables to reconstruct the dmatpawu

        matrix_i (integers) is a matrix natpawu x ndim with entries are 0 or 1
        matrix_d (deltas) is a matrix natpawu x ndim with entries are [0, 0.5)
        P (matrices) is a set of matrices natpawu x ndim x ndim
        Those three matrices allow to reconstruct the variable 'dmatpawu' used by ABINIT

        :return:
        """
        matrices_defined = []

        matrix_i = self.num_indep_matrices * [None]
        matrix_d = self.num_indep_matrices * [None]
        euler = self.num_indep_matrices * [None]

        for i in range(self.num_indep_matrices):
            nelect = self.num_electrons_dftu[i]
            val = [x for x in list(itertools.product(range(2), repeat=self.ndim)) if sum(x) == nelect]
            ii = val[np.random.randint(len(val))]
            dd = np.zeros(self.ndim)
            matrix_i[i] = list(ii)
            matrix_d[i] = list(dd)
            matrices_defined.append(self.connections[i])
            p = gram_smith_qr(self.ndim)
            euler[i] = gea_all_angles(p)

        data = {'euler_angles': euler, 'occupations': matrix_i, 'deltas': matrix_d,
                'num_matrices': self.num_indep_matrices}

        return self.new_entry(data), None

    def cross(self, ids):
        """
        Crossing algorithm used notably by GA to mix the information from several candidates
        Not implemented

        :param ids:
        :return:
        """
        pass

    def evaluate_entry(self, entry_id):
        """
        Evaluation externalized, no implemented

        :param entry_id:
        :return:
        """
        pass

    def from_dict(self, population_dict):
        pass

    def new_entry(self, properties, active=True):
        """
        Creates a new entry on the population database from given data.

        :param data: dictionary with 3 keys 'D' for deltas, 'I' for the integers
                        and eigen for the rotation matrix applied to the orbitals
        :param active: if True, the entry is enabled on the DB to be evaluated.
        :return:

        """
        status = {self.tag: active}
        entry = {'structure': self.structure.to_dict, 'properties': properties, 'status': status}
        entry_id = self.insert_entry(entry)
        pcm_log.debug('Added new entry: %s with tag=%s: %s' % (str(entry_id), self.tag, str(active)))
        return entry_id

    def is_evaluated(self, entry_id):

        """
        One candidate is considered evaluated if it contains any finite value of energy on the properties.energy field

        :param entry_id:
        :return:
        """
        entry = self.get_entry(entry_id, {'_id': 0, 'properties': 1})
        if entry['properties']['energy'] is not None:
            return True
        else:
            return False

    def check_duplicates(self, ids):
        """
        For a given list of identifiers 'ids' checks the values for the function 'distance' and return a dictionary
          where each key is the identifier of a unique candidate and the value is a list of identifiers considered
          equivalents to it.

        :param ids:  List of identifiers for wich the check will be performed
        :return:
        """
        ret = {}
        for i in range(len(ids)):
            entry_i = self.get_entry(ids[i])
            for j in range(i + 1, len(ids)):
                entry_j = self.get_entry(ids[j])
                if self.distance(ids[i], ids[j]) < 1E-3:
                    if entry_i in ret:
                        ret[entry_i].append(entry_j)
                    else:
                        ret[entry_i] = [entry_j]

    def distance(self, entry_id, entry_jd):
        """
        Measure of distance for two entries with identifiers 'entry_id' and 'entry_jd'
        TODO: The definition must be changed and compare the resulting dmatpawu instead of
        individual components

        :param entry_id: Identifier of first entry
        :param entry_jd: Identifier of second entry
        :return:
        """
        entry_i = self.get_entry(entry_id)
        entry_j = self.get_entry(entry_jd)
        dmat_i = entry_i['properties']['P']
        dmat_j = entry_j['properties']['P']
        dist_p = np.linalg.norm(dmat_j - dmat_i)
        dmat_i = entry_i['properties']['d']
        dmat_j = entry_j['properties']['d']
        dist_d = np.linalg.norm(dmat_j - dmat_i)
        return dist_d + dist_p

    def move_random(self, entry_id, factor=0.2, in_place=False, kind='move'):
        """
        Move one candidate with identifier 'entry_id' randomly with a factor
        given by 'factor'

        :param entry_id: Identifier of entry
        :param factor: Factor use to scale the randomness of change
        :param in_place: If True the candidate is changed keeping the identifier unchanged
        :param kind: Use when several algorithms are used for movement. One implemented here
        :return:
        """
        pass
        # entry_i = self.get_entry(entry_id)
        # dmat_i = entry_i['properties']['dmatpawu']
        # dmat_i += factor*np.random.random_sample(len(dmat_i))
        # self.pcdb.db.pychemia_entries.update({'_id': entry_id}, {'$set': {'properties.dmatpawu': list(dmat_i)}})

    def move(self, entry_id, entry_jd, factor=0.2, in_place=False):
        """
        Move one candidate with identifier 'entry_id' in the direction of another candidate 'entry_jd'

        :param entry_id: Identifier of first entry (Origin)
        :param entry_jd: Identifier of second entry (Target)
        :param factor: Scale factor for change, 0 scale is the 'Origin' candidate, 1 is the 'Target' candidate
                        Intermediate values will change candidates accordingly
        :param in_place: If True the candidate is changed keeping the identifier unchanged
        :return:
        """
        entry_i = self.get_entry(entry_id)
        dmat_i = np.array(entry_i['properties']['P'])
        entry_j = self.get_entry(entry_jd)
        dmat_j = np.array(entry_j['properties']['P'])
        dmat_i += factor * (dmat_j - dmat_i)

    def recover(self):
        pass

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
        entry = self.get_entry(entry_id)
        print(entry['properties']['O'], entry['properties']['D'])

    def get_duplicates(self, ids):
        return None

    def prepare_folder(self, entry_id, workdir, binary='abinit', source_dir='.'):

        if not os.path.isdir(workdir):
            os.mkdir(workdir)

        for i in ['abinit.files', 'batch.pbs']:
            if os.path.exists(workdir + os.sep + i):
                os.remove(workdir + os.sep + i)
            os.symlink(os.path.abspath(source_dir + os.sep + i), workdir + os.sep + i)

        input = InputVariables('abinit.in')
        params = self.pcdb.get_entry(entry_id)['properties']
        dmatpawu = params2dmatpawu(params, 2 * self.maxlpawu + 1)
        input['dmatpawu'] = list(dmatpawu.flatten())
        input.write(workdir + os.sep + 'abinit.in')

    def collect_data(self, entry_id, workdir):
        if os.path.isfile(workdir + '/abinit.out'):
            ao = AbinitOutput(workdir + '/abinit.out')
            if 'etot' in ao.get_energetics():
                energy = ao.get_energetics()['etot'][-1]
                print('Uploading energy data for %s' % entry_id)
                self.set_in_properties(entry_id, 'energy', energy)
                return True
            else:
                return False
        else:
            return False


def params2dmatpawu(params):
    """
    Build the variable dmatpawu from the components stored in params

    :param params: dictionary with keys 'I', 'D' and 'eigen'
    :param ndim: dimension of the correlation matrix
                5 for 'd' orbitals, 7 for 'f' orbitals
    :return:
    """
    ndim = params['ndim']

    if 'num_matrices' in params:
        num_matrices = params['num_matrices']
    else:
        num_matrices = len(params['deltas'])

    ret = np.zeros((num_matrices, ndim, ndim))

    for i in range(num_matrices):
        eigval = np.diag(params['occupations'][i]).astype(float)
        for j in range(ndim):
            if eigval[j, j] == 0:
                eigval[j, j] += params['deltas'][i][j]
            else:
                eigval[j, j] -= params['deltas'][i][j]
        rotation = gea_orthogonal_from_angles(params['euler_angles'][i])
        #print('i=%d' % i)
        #print(rotation)
        #print(eigval)
        correlation = np.dot(np.dot(rotation, eigval), rotation.T)
        ret[i] = correlation
    return ret


def dmatpawu2params(dmatpawu, ndim):
    """
    Takes the contents of the variable 'dmatpawu' and return their components as a set of 'occupations', 'deltas'
    and 'euler_angles'
    The Euler angles is a ordered list of angles that can rebuild a rotation matrix 'R'
    The rotation matrix 'R' is ensured to be an element of SO(ndim), ie det(R)=1.
    When the eigenvectors return a matrix with determinant -1 a mirror on the first dimension is applied.
    Such condition has no effect on the physical result of the correlation matrix

    :param dmatpawu: The contents of the variable 'dmatpawu'. A list of number representing N matrices ndim x ndim
    :param ndim: ndim is 5 for 'd' orbitals and 7 for 'f' orbitals
    :return:
    """
    dm = np.array(dmatpawu).reshape((-1, ndim, ndim))
    num_matrices = dm.shape[0]
    eigval = np.array([np.linalg.eigh(x)[0] for x in dm])
    occupations = np.array(np.round(eigval), dtype=int)
    deltas = np.abs(eigval - occupations)
    rotations = np.array([np.linalg.eigh(x)[1] for x in dm])

    mirror = np.eye(ndim)
    mirror[0, 0] = -1

    for i in range(len(rotations)):
        if np.linalg.det(rotations[i]) < 0:
            rotations[i] = np.dot(rotations[i], mirror)

    euler_angles = np.array([list(gea_all_angles(p)) for p in rotations])

    return {'occupations': occupations, 'deltas': deltas, 'euler_angles': euler_angles, 'num_matrices': num_matrices,
            'ndim': ndim}


def get_pattern(params, ndim):
    """

    :param params:
    :param ndim:
    :return:
    """

    eigvec = np.array(params['eigvec']).reshape((-1, ndim, ndim))
    natpawu = len(eigvec)
    connection = np.zeros((natpawu, natpawu, ndim, ndim))

    bb = np.dot(eigvec[0], np.linalg.inv(eigvec[3]))
    # connection = np.array(np.round(np.diagonal(bb)), dtype=int)

    iii = np.array(params['I'], dtype=int).reshape((-1, ndim))

    pattern = np.zeros((natpawu, natpawu))
    for i in range(natpawu):
        for j in range(i, natpawu):

            bb = np.dot(eigvec[0], np.linalg.inv(eigvec[3]))
            connection[i, j] = bb
            connection[j, i] = bb

            if np.all(np.array(iii[i] == iii[j])):
                pattern[i, j] = 1
                pattern[j, i] = 1
            else:
                pattern[i, j] = 0
                pattern[j, i] = 0

    return connection, pattern


def get_final_correlation_matrices_from_output(filename):
    rf = open(filename)
    data = rf.read()
    mainblock = re.findall('LDA\+U DATA[\s\w\d\-.=,>:]*\n\n\n', data)
    assert len(mainblock) == 1

    pattern = """For Atom\s*(\d+), occupations for correlated orbitals. lpawu =\s*([\d]+)\s*Atom\s*[\d]+\s*. Occ. for lpawu and for spin\s*\d+\s*=\s*([\d\.]+)\s*Atom\s*[\d]+\s*. Occ. for lpawu and for spin\s*\d+\s*=\s*([\d\.]+)\s*=> On atom\s*\d+\s*,  local Mag. for lpawu is[\s\d\w\.\-]*== Occupation matrix for correlated orbitals:\s*Occupation matrix for spin  1\s*([\d\.\-\s]*)Occupation matrix for spin  2\s*([\d\.\-\s]*)"""
    ans = re.findall(pattern, mainblock[0])
    print(ans)

    ret = []
    for i in ans:
        atom_data = {'atom number': int(i[0]),
                     'orbital': int(i[1]),
                     'occ spin 1': float(i[2]),
                     'occ spin 2': float(i[3])}
        matrix = [float(x) for x in i[4].split()]
        atom_data['matrix spin 1'] = list(matrix)
        matrix = [float(x) for x in i[5].split()]
        atom_data['matrix spin 2'] = list(matrix)
        ret.append(atom_data)
    return ret


def get_final_dmatpawu(filename):
    ret = get_final_correlation_matrices_from_output(filename)
    dmatpawu = []
    for i in ret:
        dmatpawu += i['matrix spin 1']
    return dmatpawu
