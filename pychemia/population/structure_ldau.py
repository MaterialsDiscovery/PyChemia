__author__ = 'Guillermo Avendano'

from _population import Population
from pychemia import pcm_log
import numpy as np
from pychemia.code.abinit import InputVariables
import os


class PopulationLDAU(Population):

    def __init__(self, name, abinit_input, natpawu, nsppol, nspinor, constrains=None):

        Population.__init__(self, name, 'global')
        if not os.path.isfile(abinit_input):
            print "Abinit input not found"
            raise ValueError
        self.input = InputVariables(abinit_input)
        self.structure = self.input.get_structure()
        self.maxlpawu = max(self.input.get_value('lpawu'))
        self.natpawu = natpawu
        self.nsppol = nsppol
        self.nspinor = nspinor
        if constrains is not None:
            assert len(constrains) == max(self.nspinor, self.nsppol)*self.natpawu
        self.constrains = constrains

    @property
    def ndim(self):
        return 2*self.maxlpawu+1

    def __str__(self):
        ret = ' Population LDA+U\n\n'
        ret += ' Name:               %s\n' % self.name
        ret += ' Tag:                %s\n' % self.tag
        ret += ' Formula:            %s\n' % self.structure.formula

        ret += ' Members:            %d\n' % len(self.members)
        ret += ' Actives:            %d\n' % len(self.actives)
        ret += ' Evaluated:          %d\n' % len(self.evaluated)
        return ret

    def add_random(self):
        """
        Creates a new set of matrices P and D such that:

        D is a diagonal matrix with entries 0 or 1
        P is a matrix whose entries are in [-1,1] and the
        sum of each column is normalized to 1

        :return:
        """
        ndim = 2*self.maxlpawu+1
        d = np.random.randint(2, size=ndim)
        P = 2*np.random.rand(ndim, ndim)-1
        c = P.sum(axis=0)
        P /= c
        return P, d

    def cross(self, ids):
        pass

    def from_dict(self, population_dict):
        pass

    def new_entry(self, data, active=True):
        properties = {'P': list(data['P'].flatten()), 'd': list(data['d'])}
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
            for j in range(i+1, len(ids)):
                entry_j = self.get_entry(ids[j])
                if self.distance(ids[i], ids[j])<1E-3:
                    if entry_i in ret:
                        ret[entry_i].append(entry_j)
                    else:
                        ret[entry_i]=[entry_j]

    def distance(self, entry_id, entry_jd):
        entry_i = self.get_entry(entry_id)
        entry_j = self.get_entry(entry_jd)
        dmat_i = entry_i['properties']['P']
        dmat_j = entry_j['properties']['P']
        dist_P = np.linalg.norm(dmat_j-dmat_i)
        dmat_i = entry_i['properties']['d']
        dmat_j = entry_j['properties']['d']
        dist_d = np.linalg.norm(dmat_j-dmat_i)
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
        dmat_i += factor*(dmat_j-dmat_i)

    def recover(self):
        pass

    def value(self, entry_id):
        pass

    def str_entry(self, entry_id):
        entry = self.get_entry(entry_id)
        print entry['properties']['P'], entry['properties']['d']

    def get_duplicates(self, ids):
        return None

def get_dmatpawu(d, P):
    return np.dot(P, np.dot(np.diag(d), np.linalg.inv(P)))

def get_P_and_d(dmatpawu):

    d, P = np.linalg.eig(dmatpawu)
    return d, P