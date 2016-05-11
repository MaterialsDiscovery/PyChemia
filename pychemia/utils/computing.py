from __future__ import print_function, unicode_literals
from builtins import str
from numbers import Number
import gzip
import hashlib
import os
import sys
import zipfile


def deep_unicode(value):
    """
    Recursively convert to unicode all string-like elements of a complex object.
    Use with care for objects that could not be converted properly to string
    This is safe for things like Atom names and symbols

    :param value: (unicode, list, dict)
    :return: (str, list, dict) The string, list or dictionary where all string-like elements
             are converted into strings
    :rtype : (str, list, dict)

    Example:
>>> deep_unicode(u'abc')
u'abc'
>>> deep_unicode([u'abc'])
[u'abc']
>>> deep_unicode({u'abc': u'def'})
{u'abc': u'def'}
>>> deep_unicode('abc')
u'abc'
    """
    if hasattr(value, 'decode'):
        value = value.decode()
    if isinstance(value, dict):
        ret = {}
        for key in value:
            ret[deep_unicode(key)] = deep_unicode(value[key])
        return ret
        # This line is Python 2.7+
        # return {unicode2string(key): unicode2string(value) for key, value in value.items()}
    elif isinstance(value, str):
        return value
    elif isinstance(value, Number):
        return value
    elif value is None:
        return None
    else:
        try:
            return [deep_unicode(element) for element in value]
        except ValueError:
            return value


def convert_color(s):
    """
    Convert a string hexadecimal representation of a color into a tuple
    RBG where each element of the tuple is between 0.0 to 1.0

    :param s: (str)
    :return: (tuple) With 3 floats representing the color in RGB
    :rtype : tuple

    Example:
>>> convert_color('FF5500')
(1.0, 0.3333333333333333, 0.0)
    """
    return float(int(s[:2], 16)) / 255, float(int(s[2:4], 16)) / 255, float(int(s[4:6], 16)) / 255


def get_int(value):
    """
    Convert a string to an integer

    :param value: (str)
    :return: (int)

    Example:
>>> get_int('3')
3
    """
    if value.isdigit():
        return int(value)
    else:
        print("ERROR: The value '%s' should be an integer" % value)
        sys.exit(2)


def get_float(value):
    """
    Convert a string to a float number

    :param value: (str)
    :return: (float)

    Example:
>>> get_float('3.0')
3.0
    """
    try:
        ret = float(value)
    except ValueError:
        raise ValueError("Could not convert '%s' into a float number" % value)
    return ret


def hashfile(filename):
    """
    Get the MD5 hash sum of a given file

    :param filename: (str)
    :return: (str)

>>> import tempfile
>>> a = tempfile.NamedTemporaryFile()
>>> hashfile(a.name)
'd41d8cd98f00b204e9800998ecf8427e'
>>> a = tempfile.NamedTemporaryFile('w')
>>> tmp= a.file.write(128000*'GAF')
>>> a.file.flush()
>>> hashfile(a.name)
'7b8a4f8a3ce222580765d577df78b782'
    """
    blocksize = 65536
    hasher = hashlib.md5()
    with open(filename, 'rb') as afile:
        buf = afile.read(blocksize)
        while len(buf) > 0:
            hasher.update(buf)
            buf = afile.read(blocksize)
    return hasher.hexdigest()


def read_file(filename):
    """
    General function to open files even if they are
    compressed

    :param filename:
    :return:

    Example:
>>> import tempfile
>>> a = tempfile.NamedTemporaryFile()
>>> read_file(a.name)
''
    """
    if not os.path.exists(filename):
        raise ValueError('Could not open file: %s' % filename)

    if filename[-3:] == '.gz':
        rf = gzip.GzipFile(filename)
    elif filename[-4:] == '.zip':
        rf = zipfile.ZipFile(filename)
    else:
        rf = open(filename)
    return rf.read()


def only_ascii(string):
    """
    Remove non-ascii characters, from monty

    :param string:
    :return:

    Example:
>>> only_ascii(u'lmnopq')
u'lmnopq'
    """
    return "".join(x for x in string if ord(x) < 128)
