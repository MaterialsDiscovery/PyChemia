from math import sqrt as _sqrt
import pychemia.utils.constants as _pc
from pychemia.utils.computing import string2number

"""
A syntactic parser for ABINIT input files ".in"
"""

__author__ = "Guillermo Avendano-Franco"
__copyright__ = "Copyright 2016"
__version__ = "1.1"
__maintainer__ = "Guillermo Avendano-Franco"
__email__ = "guillermo.avendano@uclouvain.be"
__status__ = "Development"
__date__ = "May 13, 2016"


def __parse_word(word):
    """
    Read an string and determines if it is a word or a number
    or a repetitive list os values such as 3*4.5=[4.5 , 4.5, 4.5]

    Args:
        word: An string that should be parsed

    Returns:
       result:
          Value extracted (word, number or list)

       kind:
          'word' such as 'ntime', 'eV', '*1', 'angstrom', etc
          'int' such as 1, 2, 3
          'float' such as 4.5, 6.7, etc
          'list' such as [4.5 , 4.5, 4.5]

    """
    result = None
    kind = None
    if word[0].isalpha() and word[:4] != 'sqrt' and word[:5] != '-sqrt':
        result = word
        kind = 'word'

    elif word[:4] == 'sqrt':
        result = _sqrt(float(word[5:-1]))
        kind = 'float'

    elif word[:5] == '-sqrt':
        result = -_sqrt(float(word[6:-1]))
        kind = 'float'

    elif word[0] == '*':
        result = word
        kind = 'word'

    elif word.isdigit():
        result = int(word)
        kind = "int"

    elif '*' in word:
        splt = word.split('*')

        if splt[0].isdigit():
            mult = int(splt[0])
            number, kind = string2number(splt[1])
            if number is not None:
                result = mult * [number]
                kind = 'list'
            else:
                result = None
                kind = None
    else:
        result, kind = string2number(word)

    return result, kind


def __clean_line(line):
    """
    Removes comments from a line, the comments
    could be with '#' or '!'
    Both kind of comments are used in ABINIT
    """
    line = line.strip()
    if len(line) > 0:
        if line[0] == '#' or line[0] == '!':
            result = []
        elif '#' in line:
            splt = line.split('#')
            line = splt[0]
            if '!' in line:
                splt = line.split('!')
                line = splt[0]
            result = line.split()
        elif '!' in line:
            splt = line.split('!')
            line = splt[0]
            if '#' in line:
                splt = line.split('#')
                line = splt[0]
            result = line.split()
        else:
            result = line.split()
    else:
        result = []
    return result


def __clean_input(filename):
    """
    Take an input file and return an
    string with all the comments removed

    """
    rf = open(filename, 'r')
    result = []
    for line in rf.readlines():
        ans = __clean_line(line)
        result = result + ans
    rf.close()
    return result


def parser(filename):
    """
    Take a filename and returns
    a dictionary of ABINIT input variables

    The values on each key are lists of numbers
    The values are converted to atomic units always
    that it found other variables (eV, angstrom, etc)
    Repetitions are expanded (3*1.5 => [1.5 1.5 1.5])

    Args:
       filename:
          Path for an ABINIT input file

    Returns:
       inputvars:
          Dictionary, the keys are names
          of variables and the values are
          list of numbers.

    """
    cleanlist = __clean_input(filename)

    inputvars = {}
    keyword = cleanlist[0]
    result, kind = __parse_word(keyword)
    if kind != 'word':
        raise ValueError('ERROR(parser): The first element', keyword, 'should be a keyword\n' + filename)
    varlist = []
    for word in cleanlist[1:]:
        value, kind = __parse_word(word)
        if kind == 'word':
            assert (isinstance(value, str))
            value = str(value)
            # It could be:
            # new keyword such as ntime
            # scaling factor as eV, ha
            # indefined multiplicator such as '*1'
            if value.lower()[:6] == "angstr":
                if len(varlist) == 0:
                    print('ERROR: scaling factor', value, 'before a list of values')
                else:
                    varlist = [_pc.angstrom_bohr * x for x in varlist]
            elif value.lower() == "ev":
                if len(varlist) == 0:
                    print('ERROR: scaling factor', value, 'before a list of values')
                else:
                    varlist = [_pc.eV_Ha * x for x in varlist]
            elif value.lower() == "ry":
                if len(varlist) == 0:
                    print('ERROR: scaling factor', value, 'before a list of values')
                else:
                    varlist = [0.5 * x for x in varlist]
            elif value.lower() == "t" or value.lower() == "te":
                if len(varlist) == 0:
                    print('ERROR: scaling factor', value, 'before a list of values')
                else:
                    varlist = [_pc.BField_Tesla * x for x in varlist]
            elif value.lower() == "k":
                if len(varlist) == 0:
                    print('ERROR: scaling factor', value, 'before a list of values')
                else:
                    varlist = [_pc.kb_HaK * x for x in varlist]
            elif value.lower() == "ha" or value.lower() == "hartree":
                pass
            elif value.lower() == "bohr" or value.lower() == "au":
                pass
            elif value[0] == '*':
                varlist.append(value)
            elif keyword == 'xyzfile':
                print('FOUND case XYZFILE')
                varlist.append(value)
            else:
                # New keyword
                if len(varlist) == 0:
                    print('ERROR: new keyword "%s" with no elements for previos keyword "%s"' % (value, keyword))
                else:
                    # Store the varlist in the dictionary
                    # for the previuos keyword
                    if len(varlist) == 1:
                        inputvars[keyword] = varlist[0]
                    else:
                        inputvars[keyword] = varlist
                    # New keyword
                    keyword = value
                    varlist = []
        elif kind == 'int':
            varlist.append(value)
        elif kind == 'float':
            varlist.append(value)
        elif kind == 'list':
            varlist = varlist + value
        else:
            return None

    # Final keyword
    if len(varlist) == 0:
        print('ERROR: new keyword "%s" with no elements for previos keyword "%s" ' % (value, keyword))
    else:
        # Store the varlist in the dictionary
        # for the previuos keyword
        if len(varlist) == 1:
            inputvars[keyword] = varlist[0]
        else:
            inputvars[keyword] = varlist

    return inputvars
