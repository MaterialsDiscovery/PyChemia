import os
from ..codes import CodeInput
from ...utils.computing import string2number


class ElkInput(CodeInput):

    def __init__(self, filename='elk.in'):

        CodeInput.__init__(self)
        if os.path.isfile(filename):
            self.input_file = filename
            self.read()

    def read(self):
        """
        Reads the input file and populate the variables dictionary
        """
        rf = open(self.input_file)
        lines = rf.readlines()
        rf.close()
        # All the variables stored as a python dictionary 
        ret = {}
        # Current key
        curkey = None
        index = -1
        for iline in range(len(lines)):
            # Remove leading and trailing spaces and newline jumps
            line = lines[iline].strip()
            # For debugging parser
            # print("%3d \t%s" % (iline,line))
            # Ignore empty lines
            if len(line) == 0:
                continue
            # Ignore lines that start with "!"
            elif line[0] == '!':
                continue
            elif index >= iline:
                continue
            # Special Parsing for atoms
            elif line == 'atoms':
                curkey = 'atoms'
                index = iline + 1
                line = lines[index]
                nspecies = int(line.split()[0])
                ret['atoms'] = {'order': [], 'nspecies': nspecies}
                for j in range(nspecies):
                    index += 1
                    line = lines[index].strip()
                    spfname = ''
                    for i in range(1, len(line)):
                        if line[i] == "'":
                            break
                        spfname += line[i]
                    ret['atoms']['order'].append(spfname)
                    index += 1
                    line = lines[index].strip()
                    natoms = int(line.split()[0])
                    ret['atoms'][spfname] = {'atposl': [], 'bfcmt': [], 'natoms': natoms}
                    for k in range(natoms):
                        index += 1
                        line = lines[index].strip()
                        atposl = [float(x) for x in line.split()[:3]]
                        ret['atoms'][spfname]['atposl'].append(atposl)
                        if len(line.split()) >= 6 and line.split()[4][0] != ':':
                            bfcmt = [float(x) for x in line.split()[3:6]]
                        else:
                            bfcmt = [0.0, 0.0, 0.0]
                        ret['atoms'][spfname]['bfcmt'].append(bfcmt)
            # When Line starts with an alpha character
            elif line[0].isalpha() or line == '.true.' or line == '.false.':
                # print("alph>%s" % line)
                if line == 'true' or line == '.true.':
                    ret[curkey].append(True)
                elif line == 'false' or line == '.false.':
                    ret[curkey].append(False)
                elif curkey is not None and curkey in ret:
                    if len(ret[curkey]) == 1:
                        ret[curkey] = ret[curkey][0]
                    curkey = line.split()[0]
                    ret[curkey] = []
                elif curkey is None:
                    curkey = line.split()[0]
                    ret[curkey] = []
                else:
                    print("ERROR: %s" % line)
            # When Line starts with a number
            elif line[0].isdigit() or (line[0] == '-' and line[1].isdigit()):
                # print("numb>%s" % line)
                for itoken in line.split():
                    if itoken[0] == ':':
                        break
                    elif itoken[0].isdigit() or (itoken[0] == '-' and itoken[1].isdigit()):
                        number, kind = string2number(itoken)
                        if number is None:
                            raise ValueError("Token could not be converted to number:%s" % itoken)
                        ret[curkey].append(number)
            # When line starts with single quotes
            elif line[0] == "'":
                # print("quot>%s" % line)
                word = ''
                for i in range(1, len(line)):
                    if line[i] == "'":
                        break
                    word += line[i]
                ret[curkey].append(word)
            else:
                print("Could not parse:%s" % line)
        self.variables = ret

    def __str__(self):
        """
        String representation of the input file.
        This string is used to show it to the user in interactive mode
        and to write it into a file with the method "write"

        """
        ret = ""
        for ikey in self.variables.keys():
            ret += '\n' + ikey + '\n'
            value = self.variables[ikey]
            if ikey == 'atoms':
                nspecies = value['nspecies']
                ret += " %d\n" % nspecies
                for i in range(nspecies):
                    spfname = value['order'][i]
                    ret += " '%s'\n" % spfname
                    natoms = value[spfname]['natoms']
                    ret += " %d\n" % natoms
                    for j in range(natoms):
                        ret += " %9.6f %9.6f %9.6f" % tuple(value[spfname]['atposl'][j])
                        ret += " %9.6f %9.6f %9.6f\n" % tuple(value[spfname]['bfcmt'][j])
            elif ikey == 'avec':
                for j in range(3):
                    ret += " %9.6f %9.6f %9.6f\n" % tuple(value[j * 3:j * 3 + 3])
            elif ikey == 'tasks':
                for itask in value:
                    ret += " %d\n" % itask
            elif ikey == 'wplot':
                ret += " %d %d %d\n" % tuple(value[:3])
                ret += " %9.3f %9.3f\n" % tuple(value[3:])
            else:
                if type(value) == list:
                    for i in value:
                        ret += " " + str(i)
                    ret += '\n'
                elif type(value) == int:
                    ret += " %d\n" % value
                elif type(value) == float:
                    ret += " %f\n" % value
                elif type(value) == str:
                    ret += " '%s'\n" % value
                elif type(value) == bool:
                    ret += "%s\n" % str(value).lower()
                else:
                    raise ValueError("Could not identify proper type for %s with type: %s" % (value, type(value)))

        return ret
