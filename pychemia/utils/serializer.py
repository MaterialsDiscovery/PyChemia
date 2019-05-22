
import json
from abc import ABCMeta, abstractmethod
import numpy as np
from pychemia.utils.computing import deep_unicode
from numbers import Integral, Real


class PyChemiaJsonable(object):
    """
    Abstract base class specifying how to convert objects from/to dictionaries.
    PyChemiaJsonable objects must implement a to_dict property and a from_dict static method.
    """
    __metaclass__ = ABCMeta

    @property
    @abstractmethod
    def to_dict(self):
        """
        A JSON representation of an object.
        """
        pass

    @classmethod
    def from_dict(cls, json_dict):
        """
        This implements a default from_dict method which supports all
        classes that simply try to recreate an object using the keys
        as arguments for the creation of the new object.

        :param json_dict: Recreate an object from its serialize form
        :return:
        """
        argstring = ''
        for key in json_dict:
            argstring += key + '=' + str(json_dict[key]) + ', '
        argstring = argstring[:-2]
        print(str(cls) + '(' + argstring + ')')
        return eval(str(cls) + '(' + argstring + ')')

    @property
    def to_json(self):
        """
        Returns a json string representation of the object.
        """
        return json.dumps(self)

    def save_to_file(self, filename):
        """
        Writes the json representation to a file.

        :param filename: (str) Filename for the json that will be created
        """
        with open(filename, "w") as f:
            json.dump(self, f)


def generic_serializer(value):
    """
    A generic serializer for very common values

    :param value:
    :return:
    """
    value = deep_unicode(value)

    if value is None:
        return None
    elif isinstance(value, dict):
        new_value = {}
        for i in value:
            new_value[i] = generic_serializer(value[i])
        return new_value
    elif hasattr(value, '__iter__'):
        return [generic_serializer(element) for element in value]
    elif isinstance(value, str):
        return value
    elif isinstance(value, Integral):
        return int(value)
    elif isinstance(value, Real):
        return float(value)
    elif isinstance(value, np.integer):
        return int(value)
    elif isinstance(value, np.float):
        return float(value)
    else:
        raise ValueError("Could not serialize this: %s of type: %s" % (value, type(value)))
