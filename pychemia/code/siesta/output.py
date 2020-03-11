from ..codes import CodeOutput
import os
import re
import json


class SiestaOutput(CodeOutput):

    def __init__(self, outputfile='siesta.out'):
        CodeOutput.__init__(self)
        self.outputfile = None
        self.output_values = None
        self.data = None
        if os.path.isfile(outputfile):
            self.outputfile = outputfile
        if self.is_finished:
            self.read()

    @property
    def is_finished(self):
        if self.outputfile is None:
            return False

        rf = open(self.outputfile)
        data = rf.read()
        rf.close()
        if data[-14:] == 'Job completed\n':
            return True
        else:
            return False

    def read(self):
        if not os.path.isfile(self.outputfile):
            raise ValueError("ERROR: Siesta outputfile not found: %s" % self.outputfile)
        rf = open(self.outputfile)
        self.data = rf.read()
        rf.close()

        subdata = re.findall("siesta: Final energy \(eV\):[\s\d\w\W]*\n\n", self.data)
        # print(subdata)
        if len(subdata) == 0:
            raise ValueError('No Final data could be retrieved')
        elif len(subdata) > 1:
            raise ValueError('ERROR: Wrong parsing of data')

        ret = {}
        for i in subdata[0].split('\n'):
            # Debugging parser
            # print('Line => %s' % i)
            if 'siesta:' in i:
                line = i.replace('siesta:', '').strip()
            else:
                line = i.strip()
            if 'Final energy' in line:
                master = line[:-1].strip()
            elif 'Atomic forces' in line:
                master = line[:-1].strip()
            elif 'Stress tensor' in line:
                master = line[:-1].strip()
            elif 'Cell volume' in line:
                ret['Cell volume'] = self.parse_line(line.split('=')[1])
            elif 'Pressure' in line:
                master = line[:-1].strip()
            elif '(Free)E+ p_basis*V_orbitals' in line:
                ret['(Free)E+ p_basis*V_orbitals'] = self.parse_line(line.split('=')[1])
            elif '(Free)Eharris+ p_basis*V_orbitals' in line:
                ret['(Free)Eharris+ p_basis*V_orbitals'] = self.parse_line(line.split('=')[1])
            elif 'Electric dipole (a.u.)' in line:
                ret['Electric dipole (a.u.)'] = self.parse_line(line.split('=')[1])
            elif 'Electric dipole (Debye)' in line:
                ret['Electric dipole (Debye)'] = self.parse_line(line.split('=')[1])
            elif 'Vacuum level (max, mean)' in line:
                ret['Vacuum level (max, mean)'] = self.parse_line(line.split('=')[1])
            elif 'Elapsed wall time (sec)' in line:
                ret['Elapsed wall time (sec)'] = self.parse_line(line.split('=')[1])
            elif line.strip() == '':
                continue
            elif line.strip()[:3] == '---' and line.strip()[-3:] == '---':
                continue
            elif "CPU execution times" in line:
                break
            elif master is not None:
                if '=' in line:
                    if master not in ret:
                        ret[master] = {}
                    key = line.split('=')[0].strip()
                    value = line.split('=')[1]
                    ret[master][key] = self.parse_line(value)
                else:
                    if master not in ret:
                        ret[master] = []
                    ret[master].append(self.parse_line(line))
        self.output_values = ret

    def parse_line(self, line):
        ret = []
        for i in line.split():
            try:
                value = int(i)
            except ValueError:
                try:
                    value = float(i)
                except ValueError:
                    value = i
            ret.append(value)
        if len(ret) == 1:
            ret = ret[0]
        return ret

    def show_parsed_data(self):
        if self.output_values is None:
            raise ValueError('No data has been parsed')
        else:
            print(json.dumps(self.output_values, sort_keys=True, indent=4, separators=(',', ': ')))
