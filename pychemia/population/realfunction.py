
import uuid
import numpy as np
import scipy.optimize
from ._population import Population
from pychemia.utils.mathematics import unit_vector


class RealFunction(Population):
    def __init__(self, function, ndim, limits, local_minimization=False):
        """
        Creates a simple population of points in R^N with
        N=ndim the dimensionality of space and
        using an univaluated function 'function'

        :param function: Routine to evaluate a function
        :param ndim: (int) Dimensions of function space
        :param limits: (numpy.ndarray)
        :return:
        """
        Population.__init__(self, 'Euclidean', 'global', use_mongo=False, distance_tolerance=1E-3)
        self.tag = 'global'
        self.name = 'Real Function'
        self.function = function
        self.ndim = ndim
        if len(limits) == 2:
            self.limits = np.zeros((ndim, 2))
            self.limits[:, 0] = limits[0]
            self.limits[:, 1] = limits[1]
        else:
            self.limits = np.array(limits)
            assert (self.limits.shape == (2, ndim))

        self._members = []
        self._actives = []
        self.db = {}
        self.moves = {}
        self.pcdb = None
        self.local_minimization = local_minimization

    def __str__(self):
        ret = ' Euclidean Population\n\n'
        ret += ' Name:               %s\n' % self.name
        ret += ' Tag:                %s\n' % self.tag
        ret += '\n'
        ret += ' Members:            %d\n' % len(self.members)
        ret += ' Actives:            %d\n' % len(self.actives)
        ret += ' Evaluated:          %d\n' % len(self.evaluated)
        return ret

    def new_entry(self, data, active=True):
        ident = self.new_identifier()
        x = np.atleast_1d(data)
        self.db[ident] = {'x': x, 'fx': None}
        self.evaluate_entry(ident)
        if active:
            self.actives.append(ident)
        self.members.append(ident)
        return ident

    def add_random(self):
        x = np.random.random_sample(self.ndim)
        x = x * (self.limits[:, 1] - self.limits[:, 0]) + self.limits[:, 0]
        return self.new_entry(x), None

    def coordinate(self, i):
        return self.db[i]['x']

    def cross(self, ids):
        assert len(ids) == 2
        parent1 = self.coordinate(ids[0])
        parent2 = self.coordinate(ids[1])
        if self.ndim == 1:
            son1 = 0.5 * (parent1 + parent2)
            son2 = np.abs(parent1 - parent2)
        elif self.ndim == 2:
            son1 = np.array([parent1[0], parent2[1]])
            son2 = np.array([parent2[0], parent1[1]])
        else:
            split = np.random.randint(1, self.ndim - 1)
            son1 = np.concatenate((parent1[:split], parent2[split:]))
            son2 = np.concatenate((parent2[:split], parent1[split:]))
        new_ident1 = self.new_identifier()
        self.members.append(new_ident1)
        new_ident2 = self.new_identifier()
        self.members.append(new_ident2)
        self.db[new_ident1] = {'x': son1, 'fx': None}
        self.evaluate_entry(new_ident1)
        self.db[new_ident2] = {'x': son2, 'fx': None}
        self.evaluate_entry(new_ident2)
        if self.db[new_ident1]['fx'] > self.db[new_ident2]['fx']:
            return new_ident2, new_ident1
        else:
            return new_ident1, new_ident2

    def distance(self, imember, jmember):
        # The trivial metric
        x1 = self.db[imember]['x']
        x2 = self.db[jmember]['x']
        return np.linalg.norm(x2 - x1)

    def enable(self, ident):
        if ident not in self.actives:
            self.actives.append(ident)

    def evaluate(self):
        for ident in self.members:
            self.evaluate_entry(ident)

    def evaluate_entry(self, ident):
        x = self.db[ident]['x']
        if self.local_minimization:
            localmin = scipy.optimize.minimize(self.function, x)
            if self.is_inside(localmin.x):
                self.db[ident]['x'] = localmin.x
                self.db[ident]['fx'] = localmin.fun
            else:
                self.db[ident]['fx'] = self.function(x)
        else:
            self.db[ident]['fx'] = self.function(x)

    # def evaluator_daemon(self):
    #
    #     def worker(db, function, d):
    #         while True:
    #             for entry_id in db:
    #                 if db[entry_id]['fx'] is None:
    #                     self.evaluate_entry(entry_id)
    #             for entry_id in self.db:
    #                 d[entry_id] = self.db[entry_id]
    #             time.sleep(5)
    #
    #     manager = Manager()
    #     d = manager.dict()
    #
    #     p = Process(target=worker, args=(self.db, self.function, d))
    #     p.start()
    #     return p, d

    def is_evaluated(self, i):
        if i in self.db and self.db[i]['fx'] is not None:
            return True
        else:
            return False

    def from_dict(self, population_dict):
        pass

    def disable(self, ident):
        if ident in self.actives:
            self.actives.remove(ident)

    def get_values(self, selection):
        ret = {}
        for i in selection:
            ret[i] = self.value(i)
        return ret

    def member_str(self, imember):
        ret = '('
        for i in range(self.ndim):
            ret += '%5.2f' % self.db[imember]['x'][i]
            if i < self.ndim - 1:
                ret += ', '
            else:
                ret += ') -> '
        if self.value(imember) is not None:
            ret += '%5.2f' % self.value(imember)
        else:
            ret += 'None'
        return ret

    def move(self, imember, jmember, factor=0.2, in_place=False):
        """
        Moves imember in the direction of jmember
        If in_place is True the movement occurs on the
        same address as imember

        :param factor:
        :param imember:
        :param jmember:
        :param in_place:
        :return:
        """

        x1 = self.db[imember]['x']
        x2 = self.db[jmember]['x']
        newx = x1 + factor * (x2 - x1)
        if not in_place:
            new_ident = self.new_identifier()
            self.actives.append(new_ident)
            self.members.append(new_ident)
        else:
            new_ident = imember

        # print 'Moving',imember,'located at',x1
        if new_ident not in self.moves:
            # print 'No previous'
            self.moves[new_ident] = np.vstack((x1, newx))
        else:
            # print 'With previous'
            self.moves[new_ident] = np.vstack((self.moves[new_ident], newx))
        # print self.moves[new_ident]
        self.db[new_ident] = {'x': newx, 'fx': None}
        self.evaluate_entry(new_ident)
        return new_ident

    def is_inside(self, x):
        outside = False
        for i in range(self.ndim):
            if self.limits[i, 0] > x[i] or x[i] > self.limits[i, 1]:
                outside = True
                # print('New position out of limits', x, self.limits)
        return not outside

    def move_random(self, imember, factor=0.2, in_place=False, kind='move'):

        x = np.array(self.db[imember]['x'])
        newx = x
        outside = True
        while outside:
            dx = 2 * np.random.random_sample(self.ndim) - 1
            # print 'Random movement', dx, factor
            dx = unit_vector(dx)
            newx = x + factor * dx
            outside = not self.is_inside(newx)

        if not in_place:
            new_ident = self.new_identifier()
            self.actives.append(new_ident)
            self.members.append(new_ident)
        else:
            new_ident = imember
        self.db[new_ident] = {'x': newx, 'fx': None}
        self.evaluate_entry(new_ident)
        return new_ident

    @staticmethod
    def new_identifier():
        return str(uuid.uuid4())[-12:]

    def random_population(self, n):
        for i in range(n):
            self.add_random()

    def replace_failed(self):
        pass

    def recover(self):
        pass

    def save(self):
        wf = open('population.dat', 'w')
        for i in sorted(self.members):
            wf.write("%15s %12.3f %12.3f\n" % (i, self.db[i]['x'][0], self.db[i]['x'][1]))
        wf.close()
        wf = open('members.dat', 'w')
        for i in sorted(self.members):
            wf.write("%15s\n" % i)
        wf.close()

    def save_info(self):
        pass

    def set_value(self, i, y):
        self.db[i]['fx'] = y

    def str_entry(self, entry_id):
        ret = 'x = ['
        for i in self.db[entry_id]['x']:
            ret += '%7.2e ,' % i
        ret = ret[:-1] + '] f(x) = %7.2e' % self.db[entry_id]['fx']
        return ret

    def value(self, imember):
        return self.db[imember]['fx']

    def write_change(self, change):
        # print 'Write Change', change
        if change['change'] == 'promoted':
            self.db[change['from']]['from'] = change['from']
        else:
            self.db[change['to']]['from'] = change['from']

    @property
    def actives(self):
        return self._actives

    @property
    def members(self):
        return self._members
