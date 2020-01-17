import unittest
import numpy as np
from pychemia.utils.serializer import generic_serializer


class SerializerTest(unittest.TestCase):

    def test_serializer(self):
        """
        Test (pychemia.utils.serializer)                            :
        """
        a = np.array([1, 2, 3])
        assert generic_serializer(a) == [1, 2, 3]
        b = np.array([[1, 2, 3], [4, 5, 6]])
        assert generic_serializer(b) == [[1, 2, 3], [4, 5, 6]]
        c = np.array([[[1, 2, 3], [4, 5, 6]], [[7, 8, 9], [0, 1, 2]]])
        assert generic_serializer(c) == [[[1, 2, 3], [4, 5, 6]], [[7, 8, 9], [0, 1, 2]]]
        mydict = {'a': a, 'b': b, 'c': c}
        assert generic_serializer(mydict) == {u'a': [1, 2, 3],
                                              u'b': [[1, 2, 3], [4, 5, 6]],
                                              u'c': [[[1, 2, 3], [4, 5, 6]], [[7, 8, 9], [0, 1, 2]]]}
