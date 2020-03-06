import os
from scipy.io import netcdf_file
from subprocess import Popen, PIPE


class Multibinit:
    """
    Create a Multibinit object from a given multibinit "files" file. From the files it interpret the location and read
    the contents for the several input and output files for multibinit executions.
    """
    # The default names coms from Documentation:
    # https://docs.abinit.org/guide/multibinit/
    basedir = "."
    files = "multibinit.files"
    files_list = {'in': 'multibinit.in',
                  'out': 'multibinit.out',
                  'model': 'model.XML',
                  'anharmonic': 'model_anharmonic.XML',
                  'training': 'training_set_HIST.nc'}
    input = {}

    def __init__(self, files='multibinit.files'):
        """
        Creates an Multibinit object readin the "files" file
        """
        self.files = os.path.basename(files)
        if os.path.exists(files):
            self.basedir = os.path.dirname(os.path.abspath(files))
            self.read_files()
        if os.path.isfile(self.basedir + os.sep + self.files_list['in']):
            self.read_input()

    def read_files(self):
        if not os.path.isfile(self.basedir + os.sep + self.files):
            raise ValueError("Multibinit 'files' file could not be read.")
        with open(self.basedir + os.sep + self.files) as rf:
            data = rf.readlines()
            self.files_list['in'] = data[0].strip()
            self.files_list['out'] = data[1].strip()
            self.files_list['model'] = data[2].strip()
            self.files_list['anharmonic'] = data[3].strip()
            self.files_list['training'] = data[4].strip()
            if self.files_list['anharmonic'].lower() == 'no':
                self.files_list['anharmonic'] = None
            if self.files_list['training'].lower() == 'no':
                self.files_list['training'] = None

    def read_input(self):
        self.input = {}
        rf = open(self.basedir + os.sep + self.files_list['in'])
        data = rf.readlines()
        for iline in data:
            if "#" in iline:
                iline = iline.split("#")[0].strip()
            if '=' in iline:
                variable = iline.split('=')[0].strip()
                raw_value = iline.split('=')[1].strip()
                if variable in self.input:
                    raise ValueError("Duplicate variable '%s' in input file")
                if ' ' in raw_value:
                    value_list = raw_value.split()
                    try:
                        value = [int(x) for x in value_list]
                    except ValueError:
                        try:
                            value = [float(x) for x in value_list]
                        except ValueError:
                            value = value_list
                else:
                    try:
                        value = int(raw_value)
                    except ValueError:
                        try:
                            value = float(raw_value)
                        except ValueError:
                            value = raw_value

                self.input[variable] = value

    def read_training(self):
        if self.files_list['training'] is not None:
            f = netcdf_file(self.basedir + os.sep + self.files_list['training'])

    def write_input(self, filename):
        with open(filename, 'w') as wf:
            for ivariable in sorted(self.input.keys()):
                value = self.input[ivariable]
                if isinstance(value, list):
                    wf.write("%s = " % ivariable)
                    for i in value:
                        wf.write("%s " % i)
                    wf.write("\n")
                elif isinstance(value, int):
                    wf.write("%s = %d\n" % (ivariable, value))
                elif isinstance(value, float):
                    wf.write("%s = %f\n" % (ivariable, value))
                else:
                    wf.write("%s = %s\n" % (ivariable, value))

    def write_files(self, filename):
        with open(filename, 'w') as wf:
            for ikey in ['in', 'out', 'model', 'anharmonic', 'training']:
                if ikey in ['anharmonic', 'training'] and self.files_list[ikey] is None:
                    wf.write('no\n')
                else:
                    wf.write(self.files_list[ikey] + '\n')

    def run(self, executable='multibinit', nparal=4, output_file=None, show_output=False):

        if nparal < 2:
            command = "%s < %s" % (executable, self.files)
        else:
            command = "mpirun -np %d %s < %s" % (nparal, executable, self.files)

        with Popen(command, shell=True, stdout=PIPE) as proc:
            ret = proc.stdout.read()
        if show_output:
            print(ret.decode())
        if output_file is not None:
            wf = open(output_file, 'w')
            wf.write(ret.decode())
            wf.close()
        return ret.decode()
