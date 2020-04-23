import os
import re
import time
import itertools
import subprocess
import random
import numpy as np
from .._population import Population
from pychemia import pcm_log
from pychemia.code.abinit import AbinitInput, AbinitOutput
from pychemia.utils.mathematics import gram_smith_qr, gea_all_angles, gea_orthogonal_from_angles, unit_vector
from pychemia.utils.netcdf import netcdf2dict
from pychemia.utils.serializer import generic_serializer as gs
from pychemia.runner import get_jobs, PBSRunner


class OrbitalDFTU(Population):
    def __init__(self, name, input_path='abinit.in', num_electrons_dftu=None, num_indep_matrices=None,
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

        :param input_path: The abinit input file, all the variables will be preserve for all new candidates, except
                        for 'dmatpawu' the only variable that changes.

        :param num_electrons_dftu: Example [5, 1, 5, 1]

        """
        # Call the parent class initializer to link the PyChemiaDB that will be used
        Population.__init__(self, name, 'global', distance_tolerance=0.1)
        # Checking for existence of 'abinit.in'
        if not os.path.isfile(input_path):
            raise ValueError("Abinit input not found")
        # Reading the input variables and getting the structure
        self.input_path = input_path
        self.input = AbinitInput(input_path)
        if 'dmatpawu' not in self.input.variables:
            raise ValueError("Abinit input file: %s does not contain 'dmatpawu' variable" % self.input_path)
        self.structure = self.input.get_structure()

        print('[%s] Orbital population' % self.name)
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
            # Default from ABINIT
            self.nspinor = 1

        # nspden is the number of spin-density components
        if self.input.has_variable('nspden'):
            self.nspden = self.input.get_value('nspden')
        else:
            # Default from ABINIT
            self.nspden = self.nsppol

        self.nmatrices = get_num_matrices_per_atom(self.nsppol, self.nspden, self.nspinor) * self.natpawu

        print('Variables controling the total number of matrices')
        print('nsppol : %d' % self.nsppol)
        print('nspinor: %d' % self.nspinor)
        print('nspden : %d' % self.nspden)
        print('Total number of matrices expected on dmatpawu: %d' % self.nmatrices)

        if num_electrons_dftu is None:
            abiinput = AbinitInput(input_path)
            dmatpawu = np.array(abiinput['dmatpawu']).reshape((-1, self.ndim, self.ndim))
            lpawu = abiinput['lpawu']
            maxl = max(lpawu)
            dim = 2 * maxl + 1
            params = dmatpawu2params(dmatpawu, dim)
            self.num_electrons_dftu = np.apply_along_axis(sum, 1, params['occupations'])
        else:
            self.num_electrons_dftu = np.array(num_electrons_dftu)
        print('Number of electrons for each correlation matrix: %s' % self.num_electrons_dftu)

        if num_indep_matrices is not None:
            self.num_indep_matrices = num_indep_matrices
        else:
            self.num_indep_matrices = self.nmatrices

        self.num_electrons_dftu = [int(x) for x in self.num_electrons_dftu]
        print('Number of independent matrices: %d' % self.num_indep_matrices)

        if connections is not None:
            self.connections = list(connections)
            if len(self.connections) != self.nmatrices:
                raise ValueError('Number of connections between matrices is not consistent with the number of matrices '
                                 'defined on dmatpawu')
            print('Connections: %s' % self.connections)
        else:
            self.connections = list(range(self.nmatrices))
        print("")

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

    @property
    def to_dict(self):
        ret = super(self.__class__, self).to_dict
        ret['input_path'] = self.input_path
        ret['num_electrons_dftu'] = list(self.num_electrons_dftu)
        ret['num_indep_matrices'] = self.num_indep_matrices
        ret['connections'] = list(self.connections)
        return ret

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
                'num_matrices': self.num_indep_matrices, 'ndim': self.ndim}

        return self.new_entry(data), None

    def cross(self, ids):
        """
        Crossing algorithm used notably by GA to mix the information from several candidates
        This crossing algorithm is mixing the angles of two correlation matrices preserving the
        ordering of the atoms where the angles are applied. The occupations and deltas are also
        mixed independently of the euler angles.

        :param ids:
        :return:
        """
        assert len(ids) == 2

        properties1 = self.get_correlation_params(ids[0])
        properties2 = self.get_correlation_params(ids[1])

        euler_angles1 = properties1['euler_angles']
        euler_angles2 = properties2['euler_angles']

        occupations1 = properties1['occupations']
        occupations2 = properties2['occupations']

        deltas1 = properties1['deltas']
        deltas2 = properties2['deltas']

        newdata1 = {'euler_angles': [], 'occupations': [], 'deltas': [], 'num_matrices': self.num_indep_matrices,
                    'ndim': self.ndim}
        newdata2 = {'euler_angles': [], 'occupations': [], 'deltas': [], 'num_matrices': self.num_indep_matrices,
                    'ndim': self.ndim}

        for i in range(self.num_indep_matrices):
            rnd = random.randint(0, 1)
            if rnd == 0:
                newdata1['euler_angles'].append(euler_angles1[i])
                newdata2['euler_angles'].append(euler_angles2[i])
            else:
                newdata1['euler_angles'].append(euler_angles2[i])
                newdata2['euler_angles'].append(euler_angles1[i])
            rnd = random.randint(0, 1)
            if rnd == 0:
                newdata1['occupations'].append(occupations1[i])
                newdata2['occupations'].append(occupations2[i])
                newdata1['deltas'].append(deltas1[i])
                newdata2['deltas'].append(deltas2[i])
            else:
                newdata1['occupations'].append(occupations2[i])
                newdata2['occupations'].append(occupations1[i])
                newdata1['deltas'].append(deltas2[i])
                newdata2['deltas'].append(deltas1[i])

        entry_id = self.new_entry(newdata1)
        entry_jd = self.new_entry(newdata2)
        return entry_id, entry_jd

    def collect_data(self, entry_id, workdir, filename='abinit.out'):

        if os.path.isfile(workdir + os.sep + filename):
            abo = AbinitOutput(workdir + os.sep + filename)
            dmatpawu = abo.get_dmatpawu()
            dmat = dmatpawu2params(dmatpawu, self.ndim)

            if 'etot' in abo.get_energetics():
                etot = abo.get_energetics()['etot'][-1]
                nres2 = abo.get_energetics()['nres2'][-1]
                print('Uploading energy data for %s' % entry_id)
                self.set_final_results(entry_id, dmat, etot, nres2)

                return True
            else:
                return False
        else:
            return False

    def set_final_results(self, entry_id, dmat, etot, nres2):
        self.pcdb.db.pychemia_entries.update_one({'_id': entry_id}, {'$set': {'properties.final_dmat': dmat,
                                                                          'properties.etot': etot,
                                                                          'properties.nres2': nres2}})

    def collect(self, entry_id, workdir='.'):
        idir = workdir + os.sep + str(entry_id)
        abinitout = get_final_abinit_out(idir)

        if self.get_entry(entry_id) is not None and not self.is_evaluated(entry_id):
            print("%s Setting final dmat" % entry_id)
            self.set_final_dmat(entry_id, abinitout)
        elif self.get_entry(entry_id) is None:
            print("%s Is not in database" % entry_id)
        elif self.is_evaluated(entry_id):
            print("%s Is already evaluated" % entry_id)

    def distance(self, entry_id, entry_jd):
        """
        Measure of distance for two entries with identifiers 'entry_id' and 'entry_jd'

        :param entry_id: Identifier of first entry
        :param entry_jd: Identifier of second entry
        :return:
        """

        dmat1 = self.get_correlation_params(entry_id)
        dmat2 = self.get_correlation_params(entry_jd)

        euler_angles1 = dmat1['euler_angles']
        euler_angles2 = dmat2['euler_angles']

        uvect1 = unit_vector(np.concatenate((np.cos(euler_angles1), np.sin(euler_angles1))).flatten())
        uvect2 = unit_vector(np.concatenate((np.cos(euler_angles2), np.sin(euler_angles2))).flatten())

        dist_euler = 1 - np.dot(uvect1, uvect2)

        return dist_euler

    def evaluate_entry(self, entry_id):
        """
        Evaluation externalized, no implemented

        :param entry_id:
        :return:
        """
        pass

    def evaluator(self, pbs_settings, basedir):
        print("Population size: %d" % len(self))
        while True:
            # 1. Get actives no evaluated
            ane = self.actives_no_evaluated
            if len(ane) > 0:
                print("Candidates to be evaluated: %d" % len(ane))
            else:
                print("No candidates to be evaluated")
            # 2. Get jobs in queue
            if 'user' not in pbs_settings:
                raise ValueError("PBS settings must contain a keys 'user', 'ppn' and 'walltime'")
            username = pbs_settings['user']

            jobs = get_jobs(username)
            jobnames = [jobs[x]['Job_Name'] for x in jobs]
            print("There are %d jobs in the system for user %s " % (len(jobnames), username))
            check = False
            for entry_id in ane:
                if str(entry_id) not in jobnames:
                    check = True
                else:
                    jobids = [jobs[x]['Job_Id'] for x in jobs if jobs[x]['Job_Name'] == str(entry_id)]
                    for jobid in jobids:
                        check = True
                        if jobs[jobid]['job_state'] != 'C':
                            check = False
                            break

                to_submit = False
                if check:
                    to_submit = True
                    if not os.path.isdir(basedir + os.sep + str(entry_id)):
                        self.prepare_folder(entry_id, workdir=basedir, source_dir=basedir)
                    elif os.path.isfile(basedir + os.sep + str(entry_id) + os.sep + 'COMPLETE'):
                        abinitout = get_final_abinit_out(basedir + os.sep + str(entry_id))
                        if abinitout is not None:
                            abo = AbinitOutput(abinitout)
                            # check if finished
                            if abo.is_finished:
                                # Collect results
                                self.collect(entry_id, basedir)
                                to_submit = False

                if to_submit:
                    self.submit(entry_id, basedir, pbs_settings)

            print('Sleeping for 20 minutes')
            time.sleep(1200)

    def from_dict(self, population_dict):
        return self.__init__(name=self.name,
                             input_path=population_dict['input_path'],
                             num_electrons_dftu=population_dict['num_electrons_dftu'],
                             connections=population_dict['connections'],
                             num_indep_matrices=population_dict['num_indep_matrices'])

    def get_correlation_params(self, entry_id, final=True):

        if final:
            key = 'final_dmat'
        else:
            key = 'initial_dmat'

        entry = self.get_entry(entry_id, {'properties.' + key: 1}, with_id=False)
        if entry is not None:
            ret = entry['properties'][key]
        else:
            ret = None
        return ret

    def get_final_properties(self, path):

        outputfile = get_final_abinit_out(path)
        if outputfile is None:
            raise ValueError("Not such dir: %s" % path)

        # Reading the OUTPUT
        abo = AbinitOutput(outputfile)
        dmatpawu = abo.get_dmatpawu()
        odmatpawu = np.array(dmatpawu).reshape((-1, self.ndim, self.ndim))
        oparams = dmatpawu2params(odmatpawu, self.ndim)

        if not abo.is_finished:
            print('This output is not finished')
            return None, None

        data = abo.get_energetics()
        nres2 = data['nres2'][-1]
        etot = data['etot'][-1]

        final_properties = None
        if etot is not None and oparams is not None:
            final_properties = {'etot': etot,
                                'nres2': nres2,
                                'final_dmat': gs(oparams)}

        return final_properties, etot

    def is_evaluated(self, entry_id):

        """
        One candidate is considered evaluated if it contains any finite value of energy on the properties.energy field

        :param entry_id:
        :return:
        """
        entry = self.get_entry(entry_id, {'_id': 0, 'properties': 1})
        if entry['properties']['etot'] is not None:
            return True
        else:
            return False

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
        if self.is_evaluated(entry_id):
            dmat = self.get_correlation_params(entry_id)
        else:
            dmat = self.get_correlation_params(entry_id, final=False)

        newdata = dict(dmat)
        for i in range(self.num_indep_matrices):
            for j in range(self.ndim):
                perturbation = dmat['euler_angles'][i][j] + 2 * np.random.rand() * factor - factor
                if np.pi / 2.0 > perturbation > -np.pi / 2.0:
                    newdata['euler_angles'][i][j] = perturbation

        if in_place:
            return self.update_dmat_inplace(entry_id, newdata)
        else:
            return self.new_entry(newdata, active=False)

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

        if self.is_evaluated(entry_id):
            dmat1 = self.get_correlation_params(entry_id)
        else:
            dmat1 = self.get_correlation_params(entry_id, final=False)

        if self.is_evaluated(entry_id):
            dmat2 = self.get_correlation_params(entry_jd)
        else:
            dmat2 = self.get_correlation_params(entry_jd, final=False)

        euler_angles1 = dmat1['euler_angles']
        euler_angles2 = dmat2['euler_angles']

        euler_angles_new = np.zeros((self.num_indep_matrices, int(self.ndim * (self.ndim - 1) / 2)))

        for i in range(self.num_indep_matrices):
            for j in range(int(self.ndim * (self.ndim - 1) / 2)):
                angle1 = euler_angles1[i][j]
                angle2 = euler_angles2[i][j]
                if angle1 < angle2:
                    if angle2 - angle1 < angle1 - angle2 + 2 * np.pi:
                        direction = 1  # Forward
                        angle = angle2 - angle1
                    else:
                        direction = -1  # Backward
                        angle = angle1 - angle2 + 2 * np.pi
                else:
                    if angle1 - angle2 < angle2 - angle1 + 2 * np.pi:
                        direction = -1  # Backward
                        angle = angle1 - angle2
                    else:
                        direction = 1
                        angle = angle2 - angle1 + 2 * np.pi

                euler_angles_new[i, j] = angle1 + direction * factor * angle
                if euler_angles_new[i, j] > np.pi:
                    euler_angles_new[i, j] -= -2 * np.pi
                if euler_angles_new[i, j] < -np.pi:
                    euler_angles_new[i, j] += -2 * np.pi

        newdata = dict(dmat1)
        newdata['euler_angles'] = gs(euler_angles_new)

        if in_place:
            return self.update_dmat_inplace(entry_id, newdata)
        else:
            return self.new_entry(newdata, active=False)

    def new_entry(self, dmat, active=True):
        """
        Creates a new entry on the population database from given data.

        :param dmat: density matrix params, a dictionrary with 'deltas' for deltas, 'occupations' for the integers
                        and 'euler_angles' for angles needed to recreate the rotation matrix applied to the orbitals
        :param active: if True, the entry is enabled on the DB to be evaluated.
        :return:

        """
        assert 'euler_angles' in dmat
        assert 'occupations' in dmat
        assert 'deltas' in dmat
        assert 'ndim' in dmat
        assert 'num_matrices' in dmat
        status = {self.tag: active}
        properties = {'etot': None, 'initial_dmat': dmat, 'nres2': None, 'final_dmat': None}

        entry = {'structure': self.structure.to_dict, 'properties': properties, 'status': status}
        entry_id = self.insert_entry(entry)
        pcm_log.debug('Added new entry: %s with tag=%s: %s' % (str(entry_id), self.tag, str(active)))
        return entry_id

    def plot_distance_matrix(self, filename=None, ids=None):

        import matplotlib.pyplot as plt
        if ids is None:
            ids = self.members
        if filename is None:
            filename = self.name + '.pdf'
        else:
            if filename[-3:] not in ['pdf', 'png', 'jpg']:
                raise ValueError("Filename should have extension such as pdf, png or jpg")

        m = self.distance_matrix(self.ids_sorted(ids))
        fig = plt.figure()
        ax = fig.add_axes()
        ax.imshow(m)
        ax.colorbar()
        fig.savefig(filename)
        fig.clf()

    def prepare_folder(self, entry_id, workdir='.', source_dir='.'):
        """
        Prepare directories for abinit execution

        :param entry_id:    bson.ObjectID of the entry that will be used for preparing the folder
        :param source_dir: (str) is the directory where 'abinit.files' and 'batch.pbs' should be present
                           those directories will be symbolically linked inside the individual work directories
        :param workdir: (str) Base work directory for abinit executions. Inside this directory, a set of subdirectories
                              will be created using the mongo ID as name.

        """
        # Individual workdir
        iworkdir = workdir + os.sep + str(entry_id)

        if not os.path.isdir(iworkdir):
            os.mkdir(iworkdir)

        if os.path.lexists(iworkdir + os.sep + 'abinit.files'):
            os.remove(iworkdir + os.sep + 'abinit.files')
        if not os.path.isfile(source_dir + os.sep + 'abinit.files'):
            print('WARNIG: The file %s should be present on %s, symbolic links will be created pointing '
                  'to that location' % ('abinit.files', source_dir))
        os.symlink(os.path.abspath(source_dir + os.sep + 'abinit.files'), iworkdir + os.sep + 'abinit.files')

        abiinput = AbinitInput(self.input_path)
        params = self.get_correlation_params(entry_id, final=False)
        dmatpawu = params2dmatpawu(params)
        abiinput['dmatpawu'] = list(dmatpawu.flatten())
        abiinput.write(iworkdir + os.sep + 'abinit.in')

        d_abiinput = AbinitInput(iworkdir + os.sep + 'abinit.in')
        d_dmatpawu = d_abiinput['dmatpawu']
        assert (d_dmatpawu is not None)
        d_params = dmatpawu2params(d_dmatpawu, self.ndim)
        if not np.all(np.sum(d_params['occupations'], axis=1) == np.array(self.num_electrons_dftu)):
            print('ERROR: Inconsistent number of DFT+U electrons for correlated orbitals: %s' % entry_id)
            print('From the population: %s ' % self.num_electrons_dftu)
            print('From "abinit.in":    %s ' % d_params['occupations'])

    def recover(self):
        pass

    def set_final_dmat(self, entry_id, abinitout):
        abo = AbinitOutput(abinitout)
        if not abo.is_finished:
            return None

        dmatpawu = abo.get_dmatpawu()
        if dmatpawu is None:
            return None

        odmatpawu = np.array(dmatpawu).reshape((-1, self.ndim, self.ndim))
        oparams = dmatpawu2params(odmatpawu, self.ndim)
        data = abo.get_energetics()
        if data is None:
            return None
        nres2 = data['nres2'][-1]
        etot = data['etot'][-1]

        if not np.all(np.sum(oparams['occupations'], axis=1) == np.array(self.num_electrons_dftu)):
            print('ERROR: Inconsistent number of DFT+U electrons for correlated orbitals: %s' % entry_id)
            print('From the population      : %s ' % self.num_electrons_dftu)
            print('From %20s :    %s ' % (abinitout, oparams['occupations']))

        return self.pcdb.db.pychemia_entries.update({'_id': entry_id},
                                                    {'$set': {'properties.etot': etot,
                                                              'properties.nres2': nres2,
                                                              'properties.final_dmat': gs(oparams)}})

    def str_entry(self, entry_id):
        entry = self.get_entry(entry_id, projection={'properties': 1})
        ret = "%s etot: %7.3f nres2: %9.3e" % (entry_id, entry['properties']['etot'], entry['properties']['nres2'])
        return ret

    @staticmethod
    def submit(entry_id, workdir, pbs_settings):

        if 'walltime' not in pbs_settings:
            raise ValueError('walltime is mandatory on pbs_settings')
        if 'ppn' not in pbs_settings:
            raise ValueError('ppn is mandatory on pbs_settings')
        if 'queue' not in pbs_settings:
            raise ValueError('queue is mandatory on pbs_settings')
        if 'template' not in pbs_settings:
            raise ValueError('template is mandatory on pbs_settings')
        walltime = pbs_settings['walltime']
        ppn = pbs_settings['ppn']
        queue = pbs_settings['queue']
        template = pbs_settings['template']
        if 'features' not in pbs_settings:
            features = None
        else:
            features = pbs_settings['features']
        if not os.path.isfile(template):
            raise ValueError("The file: %s must exist" % template)
        if 'pvmem' in pbs_settings:
            pvmem = pbs_settings['pvmem']
        else:
            pvmem = None
        if 'join' in pbs_settings:
            join = pbs_settings['join']
        else:
            join = None
        idir = str(entry_id)
        workdir = os.path.abspath(workdir)
        path = workdir + os.sep + idir
        print('Creating a new job at: %s' % path)
        outputs = [x for x in os.listdir(path) if x[-4:] == '.out']
        for ifile in outputs:
            os.remove(path + os.sep + ifile)
        if os.path.lexists(workdir + os.sep + 'batch.pbs'):
            os.remove(workdir + os.sep + 'batch.pbs')

        pbs = PBSRunner(workdir=path, template=template)
        pbs.set_pbs_params(nodes=1, ppn=ppn, walltime=walltime, message='ae', queue=queue, features=features,
                           join=join, pvmem=pvmem)
        pbs.write()

        jobid = pbs.submit()
        print("Entry: %s Job: %s" % (entry_id, jobid))

    def update_dmat_inplace(self, entry_id, dmat):
        return self.pcdb.db.pychemia_entries.update({'_id': entry_id},
                                                    {'$set': {'properties.initial_dmat': dmat,
                                                              'properties.etot': None,
                                                              'properties.nres2': None,
                                                              'properties.final_dmat': dmat}})

    def value(self, entry_id):
        """
        Return the energy value associated to the candidate with identifier 'entry_id'

        :param entry_id:
        :return:
        """
        entry = self.get_entry(entry_id, {'properties.etot': 1})
        return entry['properties']['etot']


def params2dmatpawu(params):
    """
    Build the variable dmatpawu from the components stored in params

    :param params: dictionary with keys 'I', 'D' and 'eigen'
    :return:
    """
    # print(list(params.keys()))

    ndim = params['ndim']

    if 'num_matrices' in params:
        num_matrices = params['num_matrices']
    else:
        num_matrices = len(params['deltas'])

    occupations = np.array(params['occupations']).reshape(num_matrices, -1)
    deltas = np.array(params['deltas']).reshape(num_matrices, -1)
    euler_angles = np.array(params['euler_angles']).reshape(num_matrices, -1)

    # print(num_matrices)
    # print(euler_angles.shape)
    # print(deltas.shape)
    # print(occupations.shape)

    ret = np.zeros((num_matrices, ndim, ndim))

    for i in range(num_matrices):
        # print(i)
        eigval = np.diag(occupations[i]).astype(float)
        for j in range(ndim):
            if eigval[j, j] == 0:
                eigval[j, j] += deltas[i, j]
            else:
                eigval[j, j] -= deltas[i, j]
        # print(eigval)
        rotation = gea_orthogonal_from_angles(euler_angles[i])
        # print(rotation)
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

    # bb = np.dot(eigvec[0], np.linalg.inv(eigvec[3]))
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


def get_num_matrices_per_atom(nsppol, nspden, nspinor):
    if nsppol == 1 and nspinor == 1 and nspden == 1:
        # Non-magnetic system (nsppol=1, nspinor=1, nspden=1):
        # One (2lpawu+1)x(2lpawu+1) dmatpawu matrix is given for each atom on which +U is applied.
        # It contains the "spin-up" occupations.
        return 1
    elif nsppol == 2 and nspinor == 1 and nspden == 2:
        # Ferromagnetic spin-polarized (collinear) system (nsppol=2, nspinor=1, nspden=2):
        # Two (2lpawu+1)x(2lpawu+1) dmatpawu matrices are given for each atom on which +U is applied.
        # They contain the "spin-up" and "spin-down" occupations.
        return 2
    elif nsppol == 1 and nspinor == 1 and nspden == 2:
        # Anti-ferromagnetic spin-polarized(collinear) system(nsppol=1, nspinor=1, nspden=2):
        # One(2lpawu + 1)x(2lpawu + 1) dmatpawu matrix is given for each atom on which +U is applied.
        # It contains the "spin-up" occupations.
        return 1
    elif nsppol == 1 and nspinor == 2 and nspden == 4:
        # Non-collinear magnetic system (nsppol=1, nspinor=2, nspden=4):
        # Two (2lpawu+1)x(2lpawu+1) dmatpawu matrices are given for each atom on which +U is applied.
        # They contains the "spin-up" and "spin-down" occupations (defined as n_up=(n+|m|)/2 and n_dn=(n-|m|)/2),
        #    where m is the integrated magnetization vector).
        return 2
    elif nsppol == 1 and nspinor == 2 and nspden == 1:
        # Non-collinear magnetic system with zero magnetization (nsppol=1, nspinor=2, nspden=1):
        # Two (2lpawu+1)x(2lpawu+1) dmatpawu matrices are given for each atom on which +U is applied.
        # They contain the "spin-up" and "spin-down" occupations;
        return 2


def get_final_abinit_out(path):
    if not os.path.isdir(path):
        raise ValueError("ERROR: Could not find folder: %s" % path)

    outputfile = None
    abos = [x for x in os.listdir(path) if x[-3:] in ['txt', 'out']]

    if len(abos) == 0:
        raise ValueError("ERROR: Not suitable ABINIT output files were found ('*out' or '*txt') for path: %s" % path)

    # Most recent mtime
    mtime = 0

    for ifile in abos:
        fpath = path + os.sep + ifile
        if os.path.getmtime(fpath) > mtime:
            abo = AbinitOutput(fpath)
            if abo.is_loaded and abo.is_finished:
                mtime = os.path.getmtime(fpath)
                outputfile = fpath

    if outputfile is None:
        print("WARNING: Could not find any good ABINIT output file: %s" % path)

    return outputfile
