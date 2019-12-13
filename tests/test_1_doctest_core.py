import doctest
import unittest
from pychemia.core.composition import Composition

from .doctest_2to3 import doctest_suite


def test_structure():
    """
    DocTests (pychemia.core.structure)                           :
    """
    import pychemia.core.structure
    dt = doctest.testmod(pychemia.core.structure, verbose=True)
    assert dt.failed == 0


def test_composition():
    """
    DocTests (pychemia.core.composition)                         :
    """
    import pychemia.core.composition
    dt = doctest.testmod(pychemia.core.composition, verbose=True)
    assert dt.failed == 0
    #suite = unittest.TestSuite()
    #suite.addTest(doctest_suite(pychemia.core.composition))
    #runner = unittest.TextTestRunner(verbosity=1)
    #result = runner.run(suite)
    #assert result.wasSuccessful()


def div(a, b):
   return a/b


class raiseTest(unittest.TestCase):
    def testraise(self):
        """
        Test exceptions in Composition                               :
        """
        a = Composition('H2O')
        self.assertRaises(ValueError, a.covalent_volume, 'torus')


if __name__ == "__main__":
    test_composition()
    test_structure()
    unittest.main()
