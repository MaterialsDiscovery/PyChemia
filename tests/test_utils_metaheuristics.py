import unittest
import pychemia.utils.metaheuristics


class MetaheuristicFunctionTest(unittest.TestCase):

    def test_metaheuristic(self):
        """
        Test (pychemia.utils.metaheuristics)                        :
        """
        for i in ['Sphere', 'Ackley', 'Rosenbrock', 'Beale', 'GoldsteinPrice', 'Booth', 'BukinN6', 'Matyas', 'LeviN13',
                  'ThreeHump', 'Easom', 'CrossInTray', 'Eggholder', 'HolderTable', 'McCormick', 'SchafferN2',
                  'SchafferN4', 'StyblinskiTang']:
            func = eval('pychemia.utils.metaheuristics.'+i+'()')

            print(i)
            print(func.mindim)
            print(func.minimum(func.mindim))
            print(func.fminimum(func.mindim))
