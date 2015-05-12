__author__ = 'Guillermo Avendano Franco'

import json
from abc import ABCMeta, abstractproperty
import numpy as np


class PyChemiaJsonable(object):
    """
    Abstract base class specifying how to convert objects from/to dictionaries.
    PyChemiaJsonable objects must implement a to_dict property and a from_dict static method.
    """
    __metaclass__ = ABCMeta

    @abstractproperty
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
        """
        argstring = ''
        for key in json_dict:
            argstring += key+'='+str(json_dict[key])+', '
        argstring = argstring[:-2]
        print str(cls)+'('+argstring+')'
        return eval(str(cls)+'('+argstring+')')

    @property
    def to_json(self):
        """
        Returns a json string representation of the object.
        """
        return json.dumps(self)

    def save_to_file(self, filename):
        """
        Writes the json representation to a file.
        """
        with open(filename, "w") as f:
            json.dump(self, f)


def generic_serializer(value):
    """
    A generic serializer for very common values

    :param value:
    :return:
    """

    if value is None:
        return None
    elif isinstance(value, np.ndarray):
        if value.ndim == 1:
            return list(value)
        elif value.ndim == 2:
            return [list(i) for i in value]
        elif value.ndim == 3:
            new_value = []
            for ivalue in value:
                new_value.append([list(i) for i in ivalue])
            return new_value
    elif isinstance(value, basestring):
        return value
    elif isinstance(value, float):
        return value
    elif isinstance(value, int):
        return value
    else:
        raise ValueError("I do not know how to covert this: ", type(value), value)
