import doctest
import unittest

from pychemia.test.doctest_2to3 import doctest_suite


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
    suite = unittest.TestSuite()
    suite.addTest(doctest_suite(pychemia.core.composition))
    runner = unittest.TextTestRunner(verbosity=1)
    result = runner.run(suite)
    assert result.wasSuccessful()


if __name__ == "__main__":
    test_composition()
