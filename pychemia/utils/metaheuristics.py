import sys
import numpy as np
from abc import ABCMeta, abstractmethod


class OptimizationTestFunction:
    __metaclass__ = ABCMeta
    """
    General class for Test Functions used for optimization
    """

    def __init__(self, mindim=1, maxdim=None, domain=np.array([-1, 1])):
        self.mindim = mindim
        self.maxdim = maxdim
        self.domain = domain

    @staticmethod
    def function(x):
        return np.sum(np.abs(x))

    @abstractmethod
    def minimum(self, ndim):
        pass

    def fminimum(self, ndim):
        x = self.minimum(ndim)
        return self.function(x)

    def get_plot_matrices(self, shape=None):
        if shape is None:
            shape = [200, 200]
        if self.domain.ndim == 1:
            dx = float(self.domain[1] - self.domain[0]) / (shape[0])
            X, Y = np.mgrid[self.domain[0]:self.domain[1]:dx, self.domain[0]:self.domain[1]:dx]
        else:
            dx = float(self.domain[0, 1] - self.domain[0, 0]) / (shape[0])
            dy = float(self.domain[1, 1] - self.domain[1, 0]) / (shape[1])
            X, Y = np.mgrid[self.domain[0, 0]:self.domain[0, 1]:dx, self.domain[1, 0]:self.domain[1, 1]:dy]
        Z = self.function(np.array([X, Y]))
        return X, Y, Z


class Sphere(OptimizationTestFunction):
    def __init__(self):
        OptimizationTestFunction.__init__(self, mindim=1, maxdim=None, domain=np.array([-5, 5]))

    @staticmethod
    def function(x):
        x = np.array(x)
        return np.sum(x.T * x.T, axis=-1).T

    def minimum(self, ndim):
        return np.zeros(ndim)


class Ackley(OptimizationTestFunction):
    def __init__(self):
        OptimizationTestFunction.__init__(self, mindim=1, maxdim=None, domain=np.array([-5, 5]))

    @staticmethod
    def function(x):
        n = len(x)
        exp1 = np.exp(-0.2 * np.sqrt(1.0 / n * np.sum(x * x)))
        exp2 = np.exp(1.0 / n * np.sum((np.cos(2 * np.pi * x)).T, axis=-1).T)
        return -20 * exp1 - exp2 + np.e + 20

    def minimum(self, ndim):
        return np.zeros(ndim)


class Rosenbrock(OptimizationTestFunction):
    def __init__(self):
        OptimizationTestFunction.__init__(self, mindim=1, maxdim=None, domain=np.array([-5, 5]))

    @staticmethod
    def function(x):
        return np.sum((100.0 * (x[1:] - x[:-1] ** 2.0) ** 2.0 + (1 - x[:-1]) ** 2.0).T, axis=-1).T

    def minimum(self, ndim):
        return np.ones(ndim)


class Beale(OptimizationTestFunction):
    def __init__(self):
        OptimizationTestFunction.__init__(self, mindim=2, maxdim=2, domain=np.array([-4.5, 4.5]))

    @staticmethod
    def function(x):
        return (1.5 - x[0] + x[0] * x[1]) ** 2 + (2.25 - x[0] + x[0] * x[1] * x[1]) ** 2 + (2.625 - x[0] + x[0] * x[1] *
                                                                                            x[1] * x[1]) ** 2

    def minimum(self, ndim):
        assert ndim == 2
        return np.array([3.0, 0.5])


class GoldsteinPrice(OptimizationTestFunction):
    def __init__(self):
        OptimizationTestFunction.__init__(self, mindim=2, maxdim=2, domain=np.array([-2, 2]))

    @staticmethod
    def function(x):
        factor1 = (19 - 14 * x[0] + 3 * x[0] ** 2 - 14 * x[1] + 6 * x[0] * x[1] + 3 * x[1] ** 2)
        factor2 = (18 - 32 * x[0] + 12 * x[0] ** 2 + 48 * x[1] - 36 * x[0] * x[1] + 27 * x[1] ** 2)
        return (1 + ((x[0] + x[1] + 1) ** 2) * factor1) * (30 + ((2 * x[0] - 3 * x[1]) ** 2) * factor2)

    def minimum(self, ndim):
        assert ndim == 2
        return np.array([0.0, -1.0])


class Booth(OptimizationTestFunction):
    def __init__(self):
        OptimizationTestFunction.__init__(self, mindim=2, maxdim=2, domain=np.array([-10, 10]))

    @staticmethod
    def function(x):
        return (x[0] + 2 * x[1] - 7) ** 2 + (2 * x[0] + x[1] - 5) ** 2

    def minimum(self, ndim):
        assert ndim == 2
        return np.array([1.0, -3.0])


class BukinN6(OptimizationTestFunction):
    def __init__(self):
        OptimizationTestFunction.__init__(self, mindim=2, maxdim=2, domain=np.array([[-15, 15], [-3, 3]]))

    @staticmethod
    def function(x):
        return 100 * np.sqrt(np.abs(x[1] - 0.01 * x[0] ** 2)) + 0.01 * np.abs(x[0] + 10)

    def minimum(self, ndim):
        assert ndim == 2
        return np.array([10.0, 1.0])


class Matyas(OptimizationTestFunction):
    def __init__(self):
        OptimizationTestFunction.__init__(self, mindim=2, maxdim=2, domain=np.array([-10, 10]))

    @staticmethod
    def function(x):
        return 0.26 * (x[0] ** 2 + x[1] ** 2) - 0.48 * x[0] * x[1]

    def minimum(self, ndim):
        assert ndim == 2
        return np.array([1.0, 1.0])


class LeviN13(OptimizationTestFunction):
    def __init__(self):
        OptimizationTestFunction.__init__(self, mindim=2, maxdim=2, domain=np.array([-10, 10]))

    @staticmethod
    def function(x):
        term1 = (np.sin(3 * np.pi * x[0])) ** 3
        term2 = (x[0] - 1) ** 2 * (1 + (np.sin(3 * np.pi * x[1])) ** 3)
        term3 = (x[1] - 1) ** 2 * (1 + (np.sin(2 * np.pi * x[1])) ** 2)
        return term1 + term2 + term3

    def minimum(self, ndim):
        assert ndim == 2
        return np.array([1.0, 1.0])


class ThreeHump(OptimizationTestFunction):
    def __init__(self):
        OptimizationTestFunction.__init__(self, mindim=2, maxdim=2, domain=np.array([-10, 10]))

    @staticmethod
    def function(x):
        return 2 * x[1] ** 2 - 1.05 * x[0] ** 4 + 1.0 / 6.0 * x[0] ** 6 + x[0] * x[1] + x[1] ** 2

    def minimum(self, ndim):
        assert ndim == 2
        return np.array([0.0, 0.0])


class Easom(OptimizationTestFunction):
    def __init__(self):
        OptimizationTestFunction.__init__(self, mindim=2, maxdim=2, domain=np.array([-100, 100]))

    @staticmethod
    def function(x):
        return -np.cos(x[0]) * np.cos(x[1]) * np.exp(-1 * ((x[0] - np.pi) ** 2 + (x[1] - np.pi) ** 2))

    def minimum(self, ndim):
        assert ndim == 2
        return np.array([np.pi, np.pi])


class CrossInTray(OptimizationTestFunction):
    def __init__(self):
        OptimizationTestFunction.__init__(self, mindim=2, maxdim=2, domain=np.array([-10, 10]))

    @staticmethod
    def function(x):
        factor1 = np.exp(np.abs(100 - np.sqrt(x[0] ** 2 + x[1] ** 2) / np.pi))
        return -1E-4 * (np.abs(np.sin(x[0]) * np.sin(x[1]) * factor1) + 1) ** 0.1

    def minimum(self, ndim):
        assert ndim == 2
        return np.array([1.34941, 1.34941])


class Eggholder(OptimizationTestFunction):
    def __init__(self):
        OptimizationTestFunction.__init__(self, mindim=2, maxdim=2, domain=np.array([-512, 512]))

    @staticmethod
    def function(x):
        return -1.0 * (x[1] + 47) * np.sin(np.sqrt(np.abs(x[1] + x[0] / 2.0 + 47))) - x[0] * np.sin(
            np.sqrt(np.abs(x[0] - x[1] - 47)))

    def minimum(self, ndim):
        assert ndim == 2
        return np.array([512, 404.2319])


class HolderTable(OptimizationTestFunction):
    def __init__(self):
        OptimizationTestFunction.__init__(self, mindim=2, maxdim=2, domain=np.array([-10, 10]))

    @staticmethod
    def function(x):
        return -1.0 * np.abs(np.sin(x[0]) * np.cos(x[1]) * np.exp(np.abs(1 - np.sqrt(x[0] ** 2 + x[1] ** 2) / np.pi)))

    def minimum(self, ndim):
        assert ndim == 2
        return np.array([8.05502, 9.664559])


class McCormick(OptimizationTestFunction):
    def __init__(self):
        OptimizationTestFunction.__init__(self, mindim=2, maxdim=2, domain=np.array([[-1.5, 4], [-3, 4]]))

    @staticmethod
    def function(x):
        return np.sin(x[0] + x[1]) + (x[0] - x[1]) ** 2 - 1.5 * x[0] + 2.5 * x[1] + 1

    def minimum(self, ndim):
        assert ndim == 2
        return np.array([8.05502, 9.66459])


class SchafferN2(OptimizationTestFunction):
    def __init__(self):
        OptimizationTestFunction.__init__(self, mindim=2, maxdim=2, domain=np.array([-100, 100]))

    @staticmethod
    def function(x):
        return 0.5 + ((np.sin(x[0] ** 2 - x[1] ** 2)) ** 2 - 0.5) / (1 + 1E-3 * (x[0] ** 2 + x[1] ** 2)) ** 2

    def minimum(self, ndim):
        assert ndim == 2
        return np.zeros(2)


class SchafferN4(OptimizationTestFunction):
    def __init__(self):
        OptimizationTestFunction.__init__(self, mindim=2, maxdim=2, domain=np.array([-100, 100]))

    @staticmethod
    def function(x):
        return 0.5 + ((np.cos(np.sin(np.abs(x[0] ** 2 - x[1] ** 2)))) ** 2 - 0.5) / (1 + 1E-3 * (
            x[0] ** 2 + x[1] ** 2)) ** 2

    def minimum(self, ndim):
        assert ndim == 2
        return np.array([0, 1.25313])


class StyblinskiTang(OptimizationTestFunction):
    def __init__(self):
        OptimizationTestFunction.__init__(self, mindim=1, maxdim=None, domain=np.array([-5, 5]))

    @staticmethod
    def function(x):
        return np.sum((x ** 4 - 16 * x ** 2 + 5 * x).T, axis=-1).T / 2.0

    def minimum(self, ndim):
        return -2.903534 * np.ones(ndim)


# class Simionescu(OptimizationTestFunction):
#
#     def __init__(self):
#         OptimizationTestFunction.__init__(self, mindim=2, maxdim=2, domain=np.array([-1.25, 1.25]))
#
#     @staticmethod
#     def function(x):
#         rt = 1
#         rs = 0.2
#         n = 8
#         return np.piecewise(x,
#                      [x[0]**2 + x[1]**2 <= (rt + rs*np.cos(n*np.arctan(x[0]/x[1])))**2,
#                       x[0]**2 + x[1]**2 > (rt + rs*np.cos(n*np.arctan(x[0]/x[1])))**2], [0.1*x[0]*x[1], 1])
#
#
#     def minimum(self, ndim):
#         assert ndim == 2
#         return -0.84852813*np.ones(ndim)


def all_tests_functions():
    current_module = sys.modules[__name__]
    f = current_module.__dict__
    return [f[x]() for x in f if hasattr(f[x], '__base__') and f[x].__base__ == OptimizationTestFunction]
