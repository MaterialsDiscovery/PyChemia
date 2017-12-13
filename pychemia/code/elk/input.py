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
        rf=open(self.input_file)
        lines=rf.readlines()
        rf.close()
        # All the variables stored as a python dictionary 
        ret={}
        # Current key
        curkey=None
        for iline in lines:
            # Remove leading and trailing spaces and newline jumps
            line=iline.strip()
            # Ignore empty lines
            if len(line)==0:
                continue
            # Ignore lines that start with "!"
            elif line[0] == '!':
                continue
            # When Line starts with an alpha character
            elif line[0].isalpha():
                # print("alph>%s" % line)
                if line=='true':
                    ret[curkey].append(True)
                elif line=='false':
                    ret[curkey].append(False)
                elif curkey is not None and curkey in ret:
                    if len(ret[curkey])==1:
                        ret[curkey]=ret[curkey][0]
                    curkey=line.split()[0]
                    ret[curkey]=[]
                elif curkey is None:
                    curkey=line.split()[0]
                    ret[curkey]=[]
                else:
                    print("ERROR: %s" % line)
            # When Line starts with a number
            elif line[0].isdigit() or (line[0]=='-' and line[1].isdigit()):
                # print("numb>%s" % line)
                for itoken in line.split():
                    if itoken[0] == ':':
                        break
                    elif itoken[0].isdigit() or (itoken[0]=='-' and itoken[1].isdigit()):
                        number, kind = string2number(itoken)
                        if number is None:
                            raise ValueError("Token could not be converted to number:%s" % itoken)
                        ret[curkey].append(number)
            # When line starts with single quotes
            elif line[0]=="'":
                # print("quot>%s" % line)
                word=''
                for i in range(1, len(line)):
                    if line[i]=="'":
                        break
                    word+=line[i]
                ret[curkey].append(word)
            else:
                print("Could not parse:%s" % line)

        self.variables=ret


    def __str__(self):
        """
        String representation of the input file.
        This string is used to show it to the user in interactive mode
        and to write it into a file with the method "write"

        """
        ret = ""
        for ikey in self.variables.keys():
            ret+='\n'+ikey+'\n'
            value = self.variables[ikey]
            if ikey == 'atoms':
                ret+=" %d\n" % value[0]
                ret+=" '%s'\n" % value[1]
                natoms = value[2]
                ret+=" %d\n" % value[2]
                for j in range(natoms):
                    ret+=" %9.6f %9.6f %9.6f %9.6f %9.6f %9.6f\n" % tuple(value[3+j*6:3+j*6+6])
            elif ikey == 'avec':
                for j in range(3):
                    ret+=" %9.6f %9.6f %9.6f\n" % tuple(value[j*3:j*3+3])
            else:
                if type(value) == list:
                    for i in value:
                        ret+=" "+str(i)
                    ret+='\n'
                if type(value) == int:
                    ret+=" %d\n" % value
                if type(value) == float:
                    ret+=" %f\n" % value
                if type(value) == str:
                    ret+=" '%s'\n" % value

        return ret
