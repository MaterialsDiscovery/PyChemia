import doctest
import unittest


def broken_function():
    raise Exception('This is broken')


class MyTestCase(unittest.TestCase):
    def test(self):
        """
        Tests from exceptions on utils      :
        """
        from pychemia.utils.periodic import atomic_number
        with self.assertRaises(Exception) as context:
            atomic_number(['H', u'A'])
        self.assertTrue('Atomic symbol not found' in context.exception)

        from pychemia.utils.computing import read_file
        with self.assertRaises(Exception) as context:
            read_file('/dev/abc')
        self.assertTrue('Could not open file: /dev/abc' in context.exception)

        from pychemia.utils.computing import get_float
        with self.assertRaises(Exception) as context:
            get_float('3i')
        self.assertTrue("Could not convert '3i' into a float number" in context.exception)


def test_periodic():
    """
    Tests from doctests for periodic    :
    """
    import pychemia.utils.periodic
    dt = doctest.testmod(pychemia.utils.periodic, verbose=True)
    assert dt.failed == 0


def test_mathematics():
    """
    Tests from doctests for mathematics :
    """
    import pychemia.utils.mathematics
    dt = doctest.testmod(pychemia.utils.mathematics, verbose=True)
    assert dt.failed == 0


def test_computing():
    """
    Tests from doctests for computing   :
    """
    import pychemia.utils.computing
    dt = doctest.testmod(pychemia.utils.computing, verbose=True)
    assert dt.failed == 0


if __name__ == '__main__':
    unittest.main()
