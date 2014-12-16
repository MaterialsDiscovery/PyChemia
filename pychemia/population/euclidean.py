__author__ = 'Guillermo Avendano-Franco'

import uuid
import numpy as np


class EuclideanPopulation():

    def __init__(self, function, ndim, limits, delta=0.1):
        self.function = function
        self.ndim = ndim
        self.delta = delta
        if len(limits) == 2:
            self.limits = np.zeros((ndim, 2))
            self.limits[:, 0] = limits[0]
            self.limits[:, 1] = limits[1]
        else:
            self.limits = np.array(limits)
            assert (self.limits.shape == (2, ndim))

        self.members = []
        self.actives = []
        self.evaluated = []
        self.db = {}

    @property
    def all_entries(self):
        return self.members

    def is_evaluated(self, i):
        if i in self.db and self.db[i]['fx'] is not None:
            return True
        else:
            return False

    def coordinate(self, i):
        return self.db[i]['x']

    def set_value(self, i, y):
        self.db[i]['fx'] = y

    def random_population(self, n):
        for i in range(n):
            self.add_random()

    def check_duplicates(self):
        ret = []
        ids = self.actives
        values = np.array([self.value(i) for i in self.actives if i in self.evaluated])
        if len(values) == 0:
            return ret
        #else:
        #    print 'Values', values
        argsort = np.argsort(values)
        diffs = np.ediff1d(values[argsort])
        #print 'Values     ', values[argsort]
        #print 'Differences', diffs
        for i in range(len(values)-1):
            idiff = diffs[i]
            if idiff < 1E-2:
                ident1 = ids[argsort[i]]
                ident2 = ids[argsort[i+1]]
                distance = self.distance(ident1, ident2)
                if distance < 1E-2:
                    #print x1, x2, "are too close"
                    fx1 = self.value(ident1)
                    fx2 = self.value(ident2)
                    if fx2 < fx1:
                        ret.append(ident1)
                    else:
                        ret.append(ident2)
        return ret

    def distance(self, imember, jmember):
        # The trivial metric
        x1 = self.db[imember]['x']
        x2 = self.db[jmember]['x']
        return np.linalg.norm(x2-x1)

    @staticmethod
    def new_identifier():
        return str(uuid.uuid4())[-12:]

    def add_random(self):
        ident = self.new_identifier()
        x = np.random.rand(self.ndim)
        x = x*(self.limits[:, 1]-self.limits[:, 0])+self.limits[:, 0]
        self.db[ident] = {'x': x, 'fx': None}
        self.actives.append(ident)
        self.members.append(ident)
        return ident

    def add_modified(self, ident):
        new_ident = self.new_identifier()
        while True:
            x = self.db[ident]['x'].copy()
            for i in range(self.ndim):
                x[i] += 2*self.delta*np.random.rand() - self.delta
            if self.limits[i, 0] < x[i] < self.limits[i, 1]:
                break

        self.db[new_ident] = {'x': x, 'fx': None}
        self.actives.append(new_ident)
        self.members.append(new_ident)
        return new_ident

    def disable(self, ident):
        if ident not in self.actives:
            raise ValueError(ident + ' not in actives')
        self.actives.remove(ident)

    @property
    def fraction_evaluated(self):
        ret = np.sum([1 for i in self.actives if i in self.evaluated])
        return float(ret)/len(self.actives)

    def active_no_evaluated(self):
        ret = []
        for i in self.actives:
            if i not in self.evaluated:
                ret.append(i)
        return ret

    def value(self, imember):
        return self.db[imember]['fx']

    def save(self):
        wf = open('population.dat', 'w')
        for i in sorted(self.members):
            wf.write("%15s %12.3f %12.3f\n" % (i, self.db[i]['x'][0], self.db[i]['x'][1]))
        wf.close()
        wf = open('members.dat', 'w')
        for i in sorted(self.members):
            wf.write("%15s\n" % i)
        wf.close()

    def member_str(self, imember):
        ret = '('
        for i in range(self.ndim):
            ret += '%5.2f' % self.db[imember]['x'][i]
            if i < self.ndim-1:
                ret += ', '
            else:
                ret += ') -> '
        if self.value(imember) is not None:
            ret += '%5.2f' % self.value(imember)
        else:
            ret += 'None'
        return ret

    def move(self, imember, jmember, in_place=False):
        """
        Moves imember in the direction of jmember
        If in_place is True the movement occurs on the
        same address as imember

        :param imember:
        :param jmember:
        :param in_place:
        :return:
        """

        x1 = self.db[imember]['x']
        x2 = self.db[jmember]['x']
        uvector = (x2-x1)/np.linalg.norm(x2 - x1)
        if not in_place:
            new_ident = self.new_identifier()
            self.actives.append(new_ident)
            self.members.append(new_ident)
        else:
            new_ident = imember
        self.db[new_ident] = {'x': x1 + self.delta*uvector, 'fx': None}
        return new_ident
