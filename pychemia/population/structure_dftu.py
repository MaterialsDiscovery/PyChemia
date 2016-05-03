from __future__ import print_function
import os
import itertools
import numpy as np
from pychemia.code.abinit import InputVariables
from ._population import Population
from pychemia import pcm_log
from pychemia.utils.mathematics import gram_smith_qr


class PopulationDFTU(Population):
    def __init__(self, name, abinit_input='abinit.in', natpawu=None, oxidations=None):

        Population.__init__(self, name, 'global')
        if not os.path.isfile(abinit_input):
            print("Abinit input not found")
            raise ValueError
        self.input = InputVariables(abinit_input)
        self.structure = self.input.get_structure()
        self.maxlpawu = max(self.input.get_value('lpawu'))

        if natpawu is None:
            spinat = self.input.get_value('spinat')
            if spinat is None:
                print('I could not determine the number of atoms playing a role in the DFT+U calculation')
                raise ValueError('Could not determine natpawu')
            else:
                spinat = np.array(spinat).reshape((-1, 3))
                self.natpawu = np.sum(np.apply_along_axis(np.linalg.norm, 1, spinat) != 0)
        else:
            self.natpawu = natpawu
        if self.input.has_variable('nsppol'):
            self.nsppol = self.input.get_value('nsppol')
        else:
            self.nsppol = 1
        if self.input.has_variable('nspinor'):
            self.nspinor = self.input.get_value('nspinor')
        else:
            self.nspinor = 1

        self.oxidations = oxidations
        self.connection = [1, -1, -1, 1]

    @property
    def ndim(self):
        return 2 * self.maxlpawu + 1

    def __str__(self):
        ret = ' Population LDA+U\n\n'
        ret += ' Name:               %s\n' % self.name
        ret += ' Tag:                %s\n' % self.tag
        ret += ' Formula:            %s\n' % self.structure.formula
        ret += ' natpawu:            %d\n' % self.natpawu
        ret += ' connection:         %s\n' % self.connection

        ret += ' Members:            %d\n' % len(self.members)
        ret += ' Actives:            %d\n' % len(self.actives)
        ret += ' Evaluated:          %d\n' % len(self.evaluated)
        return ret

    def add_random(self):
        """
        Creates a new set of variables to reconstruct the dmatpawu

        I (integers) is a matrix natpawu x ndim with entries are 0 or 1
        D (deltas) is a matrix natpawu x ndim with entries are [0, 0.5)
        P (matrices) is a set of matrices natpawu x ndim x ndim
        Those three matrices allow to reconstruct the variable 'dmatpawu' used by ABINIT

        :return:
        """
        I = np.zeros((self.natpawu, self.ndim), dtype=int)
        D = np.zeros((self.natpawu, self.ndim))
        eigvec = np.zeros((self.natpawu, self.ndim, self.ndim))

        # I and D for atoms 0 and 3
        if self.oxidations[0] < 0:
            nelect = self.ndim + self.oxidations[0]
        else:
            nelect = self.oxidations[0]
        val = [x for x in list(itertools.product(range(2), repeat=self.ndim)) if sum(x) == nelect]
        ii = val[np.random.randint(len(val))]
        dd = 0.0 * np.random.random_sample((self.ndim,))
        I[0] = ii
        D[0] = dd
        I[3] = ii
        D[3] = dd

        # I and D for atoms 1 and 2
        if self.oxidations[1] < 0:
            nelect = self.ndim + self.oxidations[1]
        else:
            nelect = self.oxidations[1]
        val = [x for x in list(itertools.product(range(2), repeat=self.ndim)) if sum(x) == nelect]
        ii = val[np.random.randint(len(val))]
        dd = 0.0 * np.random.random_sample((self.ndim,))
        I[1] = ii
        D[1] = dd
        I[2] = ii
        D[2] = dd

        p = gram_smith_qr(self.ndim)
        eigvec[0] = p
        eigvec[3] = np.dot(np.diag(self.connection), p)

        p = gram_smith_qr(self.ndim)
        eigvec[1] = p
        eigvec[2] = np.dot(np.diag(self.connection), p)

        data = {'eigvec': eigvec, 'I': I, 'D': D}

        return self.new_entry(data)

    def cross(self, ids):
        pass

    def from_dict(self, population_dict):
        pass

    def new_entry(self, data, active=True):
        properties = {'eigvec': list(data['eigvec'].flatten()),
                      'D': list(data['D'].flatten()),
                      'I': list(data['I'].flatten())}
        status = {self.tag: active}
        entry_id = self.pcdb.insert(structure=self.structure, properties=properties, status=status)
        pcm_log.debug('Added new entry: %s with tag=%s: %s' % (str(entry_id), self.tag, str(active)))
        return entry_id

    def is_evaluated(self, entry_id):
        pass

    def check_duplicates(self, ids):
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
        entry_i = self.get_entry(entry_id)
        entry_j = self.get_entry(entry_jd)
        dmat_i = entry_i['properties']['P']
        dmat_j = entry_j['properties']['P']
        dist_P = np.linalg.norm(dmat_j - dmat_i)
        dmat_i = entry_i['properties']['d']
        dmat_j = entry_j['properties']['d']
        dist_d = np.linalg.norm(dmat_j - dmat_i)
        return dist_d + dist_P

    def move_random(self, entry_id, factor=0.2, in_place=False, kind='move'):
        pass
        # entry_i = self.get_entry(entry_id)
        # dmat_i = entry_i['properties']['dmatpawu']
        # dmat_i += factor*np.random.random_sample(len(dmat_i))
        # self.pcdb.db.pychemia_entries.update({'_id': entry_id}, {'$set': {'properties.dmatpawu': list(dmat_i)}})

    def move(self, entry_id, entry_jd, factor=0.2, in_place=False):
        entry_i = self.get_entry(entry_id)
        dmat_i = np.array(entry_i['properties']['P'])
        entry_j = self.get_entry(entry_jd)
        dmat_j = np.array(entry_j['properties']['P'])
        dmat_i += factor * (dmat_j - dmat_i)

    def recover(self):
        pass

    def value(self, entry_id):
        pass

    def str_entry(self, entry_id):
        entry = self.get_entry(entry_id)
        print(entry['properties']['P'], entry['properties']['d'])

    def get_duplicates(self, ids):
        return None


def params_reshaped(params, ndim):
    iii = np.array(params['I'], dtype=int).reshape((-1, ndim))
    ddd = np.array(params['D']).reshape((-1, ndim))
    eigvec = np.array(params['eigvec']).reshape((-1, ndim, ndim))
    return iii, ddd, eigvec


def params2dmatpawu(params, ndim):
    iii, ddd, eigvec = params_reshaped(params, ndim)

    eigval = np.array(iii, dtype=float)
    for i in range(len(eigval)):
        for j in range(ndim):
            if iii[i, j] == 0:
                eigval[i, j] += ddd[i, j]
            else:
                eigval[i, j] -= ddd[i, j]
    dm = np.zeros((len(eigvec), ndim, ndim))
    for i in range(len(eigvec)):
        dm[i] = np.dot(eigvec[i], np.dot(np.diag(eigval[i]), np.linalg.inv(eigvec[i])))
    return dm


def dmatpawu2params(dmatpawu, ndim):
    dm = np.array(dmatpawu).reshape((-1, ndim, ndim))
    eigval = np.array([np.linalg.eigh(x)[0] for x in dm])
    iii = np.array(np.round(eigval), dtype=int)
    ddd = np.abs(eigval - iii)
    eigvec = np.array([np.linalg.eigh(x)[1] for x in dm])

    params = {'I': list(iii.flatten()),
              'D': list(ddd.flatten()),
              'eigvec': list(eigvec.flatten())}
    return params


def get_pattern(params, ndim):
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

            if np.all(iii[i] == iii[j]):
                pattern[i, j] = 1
                pattern[j, i] = 1
            else:
                pattern[i, j] = 0
                pattern[j, i] = 0

    return connection, pattern
