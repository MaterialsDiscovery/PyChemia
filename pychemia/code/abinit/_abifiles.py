import os as _os
import subprocess as _subprocess
from ._utils import netcdf2dict, psp_name


class AbiFiles:
    """
    Read an 'abinit.files' file and extract the
    names of input, output filenames and put them
    into dictionary called 'files'.
    """
    basedir = "."
    filename = "abinit.files"
    files = {'in': 'abinit.in',
             'out': 'abinit.out',
             'tmpin': 'abinit-i',
             'tmpout': 'abinit-o',
             'tmp': 'abinit',
             'psps': ''}

    def __init__(self, *args, **kwargs):
        """
        Creates an abifiles object

        Args:
            args:
                If args[0] is an existing filename, assume that as the
                name of the 'abinit.files' and their path is the
                basedir, otherwise will consider that as the basedir.
            kwargs:
                Valid keywords are: basedir,filename,files,in.out,
                tmpin,tmpout,tmp and psps
        """
        self.inp = None
        if len(args) == 1:
            if _os.path.isfile(args[0]):
                (self.basedir, self.filename) = _os.path.split(args[0])
                if self.basedir == "":
                    self.basedir = "."
                inputfile = open(args[0], "r")
                self.files['in'] = inputfile.readline()[:-1]
                self.files['out'] = inputfile.readline()[:-1]
                self.files['tmpin'] = inputfile.readline()[:-1]
                self.files['tmpout'] = inputfile.readline()[:-1]
                self.files['tmp'] = inputfile.readline()[:-1]
                self.files['psps'] = map(str.strip, inputfile.readlines())

            elif _os.path.isdir(args[0]):
                self.basedir = args[0]

            elif not _os.path.exists(args[0]):
                self.basedir = args[0]

        if len(kwargs) > 0:
            if 'basedir' in kwargs:
                self.basedir = kwargs['basedir']
                print(self.filename)
            if 'filename' in kwargs:
                self.filename = kwargs['filename']
            if 'files' in kwargs:
                self.files['in'] = kwargs['files'] + ".in"
                self.files['out'] = kwargs['files'] + ".out"
                self.files['tmpin'] = kwargs['files'] + "-i"
                self.files['tmpout'] = kwargs['files'] + "-o"
                self.files['tmp'] = kwargs['files']
            for x in ['in', 'out', 'tmpin', 'tmpout', 'tmp', 'psps']:
                if x in kwargs:
                    self.files[x] = kwargs[x]

    def write(self, filename):
        """
        Write the file 'filename' with the format of an
        usual '.files'

        :param filename: (str) Filename to write the 'abinit.files' file
        """
        wf = open(filename, 'w')
        wf.write(self.__str__())
        wf.close()

    def create(self):
        """
        Create the '.files' file and the directory
        if it does not exists
        """

        if not _os.path.exists(self.basedir):
            _os.makedirs(self.basedir)
        self.write(self.basedir + "/" + self.filename)
        # Write the input file
        if self.inp is not None:
            self.inp.write(self.get_input_filename())

    def __str__(self):
        """
        String version of the object, shows the same
        info as the contents of the .files file
        """
        ret = ''
        for x in ['in', 'out', 'tmpin', 'tmpout', 'tmp']:
            ret = ret + self.files[x] + "\n"
        for x in self.files['psps']:
            ret = ret + x + "\n"
        return ret

    def __repr__(self):
        """
        Representation of an abifiles object
        """
        ret = "basedir:           " + self.basedir + '\n' + \
              "filename:          " + self.filename + '\n' + \
              "        \-> in     " + self.files['in'] + '\n' + \
              "        \-> out    " + self.files['out'] + '\n' + \
              "        \-> tmpin  " + self.files['tmpin'] + '\n' + \
              "        \-> tmpout " + self.files['tmpout'] + '\n' + \
              "        \-> tmp    " + self.files['tmp'] + '\n'
        for x in self.files['psps']:
            ret = ret + "        \-> psps   " + x + '\n'
        return ret

    def get_input_filename(self):
        """
        Return the input filename
        """
        return self.basedir + "/" + self.files['in']

    def get_output_filename(self):
        """
        Return the output filename
        """
        return self.basedir + "/" + self.files['out']

    def get_out_filename(self):
        """
        Return the OUT.nc filename
        """
        return self.basedir + "/" + self.files['tmpout'] + '_OUT.nc'

    def clean(self):
        """
        Remove all the output
        """
        if _os.path.isdir(self.basedir):
            _os.remove(self.get_out_filename())
            outfile = self.files['out']
            outs = [x for x in _os.listdir(self.basedir) if x[:len(outfile)] == outfile]
            for i in outs:
                _os.remove(self.basedir + '/' + i)

    def cleanall(self):
        if _os.path.isdir(self.basedir):
            outfile = self.files['out']
            outs = [x for x in _os.listdir(self.basedir)]
            for i in outs:
                _os.remove(self.basedir + '/' + i)
            self.create()

    def set_input(self, inp):
        """
        Associate an inputvars object to the abifiles object
        """
        self.inp = inp

    def get_output(self):
        """
        Return the output as a dictionary
        """
        return netcdf2dict(self.get_out_filename())

    def execute(self, abinit_binary):
        """
        Utility that copy a given script and execute the given
        command inside the directory
        """
        cwd = _os.getcwd()
        _os.chdir(self.basedir)
        abifile = open(self.filename)
        logfile = open('abinit.log', 'w')
        _subprocess.call([abinit_binary], stdin=abifile, stdout=logfile)
        logfile.close()
        abifile.close()
        _os.chdir(cwd)

    def set_psps(self, exchange='LDA', kind='FHI'):
        """
        Set the pseudopotentials acording to
        given exchange and kind
        The pair (exchange,kind) could be:
        ('LDA','FHI')
        ('LDA','TM')
        ('GGA','FHI')
        """
        if self.inp is None:
            print('ABINIT input file not declared, the pseudopotentials cannot be set')
        else:
            self.files['psps'] = []
            pspdir = _os.getenv('HOME') + '/.abinit/' + exchange + '_' + kind
            if isinstance(self.inp.variables['znucl'], (int, float)):
                lstznucl = [self.inp.variables['znucl']]
            else:
                lstznucl = self.inp.variables['znucl']

            for i in lstznucl:
                self.files['psps'].append(pspdir + '/' + psp_name(i, exchange, kind))
