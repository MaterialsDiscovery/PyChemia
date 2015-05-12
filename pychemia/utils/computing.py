__author__ = 'Guillermo Avendano-Franco'


def unicode2string(value):
    """
    Recursively convert the unicode elements of a python object into
    python strings.
    Use with care for objects that could not be converted properly to string

    :param value: (unicode, list, dict)
    :return: (str, list, dict)

    Examples

>>> unicode2string(u'abc')
'abc'
>>> unicode2string([u'abc'])
['abc']
>>> unicode2string({u'abc': u'def'})
{'abc': 'def'}
    """
    if isinstance(value, dict):
        return {unicode2string(key): unicode2string(value) for key, value in value.iteritems()}
    elif isinstance(value, list):
        return [unicode2string(element) for element in value]
    elif isinstance(value, unicode):
        return value.encode('utf-8')
    else:
        return value


def convert_color(s):
    return float(int(s[:2], 16))/255, float(int(s[2:4], 16))/255, float(int(s[4:6], 16))/255
