import doctest
import unittest
from .doctest_2to3 import doctest_suite


def broken_function():
    raise Exception('This is broken')


class MyTestCase(unittest.TestCase):
    def test(self):
        """
        DocTests (pychemia.utils) [exceptions]                       :
        """
        from pychemia.utils.periodic import atomic_number
        with self.assertRaises(Exception) as context:
            atomic_number(['H', u'A'])
        # self.assertTrue(u'Atomic symbol not found' == context.exception)

        from pychemia.utils.computing import read_file
        with self.assertRaises(Exception) as context:
            read_file('/dev/abc')
        # self.assertTrue('Could not open file: /dev/abc' in context.exception)

        from pychemia.utils.computing import get_float
        with self.assertRaises(Exception) as context:
            get_float('3i')
            # self.assertTrue("Could not convert '3i' into a float number" in context.exception)


def test_periodic():
    """
    DocTests (pychemia.utils.periodic)                           :
    """
    import pychemia.utils.periodic
    dt = doctest.testmod(pychemia.utils.periodic, verbose=True)
    assert dt.failed == 0


def test_mathematics():
    """
    DocTests (pychemia.utils.mathematics)                        :
    """
    import pychemia.utils.mathematics
    dt = doctest.testmod(pychemia.utils.mathematics, verbose=True)
    assert dt.failed == 0


def test_computing():
    """
    DocTests (pychemia.utils.computing)                          :
    """
    import pychemia.utils.computing
    suite = unittest.TestSuite()
    suite.addTest(doctest_suite(pychemia.utils.computing))
    runner = unittest.TextTestRunner(verbosity=1)
    result = runner.run(suite)
    assert result.wasSuccessful()


if __name__ == "__main__":
    unittest.main(defaultTest='test_computing')
    unittest.main()
