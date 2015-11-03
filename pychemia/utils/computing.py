import sys
import hashlib


def unicode2string(value):
    """
    Recursively convert the unicode elements of a python object into
    python strings.
    Use with care for objects that could not be converted properly to string
    This is safe for things like Atom names and symbols

    :param value: (unicode, list, dict)
    :return: (str, list, dict) The string, list or dictionary where all string-like elements
             are converted into strings
    :rtype : (str, list, dict)

    Example:
>>> unicode2string(u'abc')
'abc'
>>> unicode2string([u'abc'])
['abc']
>>> unicode2string({u'abc': u'def'})
{'abc': 'def'}
    """
    if isinstance(value, dict):
        ret = {}
        for key in value:
            ret[unicode2string(key)] = unicode2string(value[key])
        return ret
        # This line is Python 2.7+
        # return {unicode2string(key): unicode2string(value) for key, value in value.items()}
    elif isinstance(value, list):
        return [unicode2string(element) for element in value]
    elif isinstance(value, str):
        return value.encode('utf-8')
    else:
        return value


def convert_color(s):
    """
    Convert a string hexadecimal representation of a color into a tuple
    RBG where each element of the tuple is between 0.0 to 1.0

    :param s: (str)
    :return: (tuple) With 3 floats representing the color in RGB
    :rtype : tuple

>>> import pychemia
>>> pychemia.utils.computing.convert_color('FF5500')
(1.0, 0.3333333333333333, 0.0)
    """
    return float(int(s[:2], 16))/255, float(int(s[2:4], 16))/255, float(int(s[4:6], 16))/255


def get_int(value):
    if value.isdigit():
        return int(value)
    else:
        print "ERROR: The value '%s' should be an integer" % value
        sys.exit(2)


def get_float(value):
    try:
        ret = float(value)
    except ValueError:
        print "ERROR: The value '%s' should be an float number" % value
        sys.exit(2)
    return ret


def hashfile(filename):
    blocksize = 65536
    hasher = hashlib.md5()
    with open(filename, 'rb') as afile:
        buf = afile.read(blocksize)
        while len(buf) > 0:
            hasher.update(buf)
            buf = afile.read(blocksize)
    return hasher.hexdigest()