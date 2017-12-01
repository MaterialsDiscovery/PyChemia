from abc import ABCMeta, abstractmethod, abstractproperty
import os
import subprocess


class Codes:
    __metaclass__ = ABCMeta

    def __init__(self):
        self.stdin_file = None
        self.stdout_file = None
        self.stderr_file = None
        self.stdin_filename = None
        self.stdout_filename = None
        self.stderr_filename = None
        self.binary = None
        self.runner = None
        self.workdir = None

    @abstractmethod
    def initialize(self, structure, workdir=None, kpoints=None, binary=None):
        pass

    @abstractmethod
    def set_inputs(self):
        pass

    @abstractmethod
    def get_outputs(self):
        pass

    @abstractmethod
    def finalize(self):
        pass

    def run(self, use_mpi=False, omp_max_threads=0, mpi_num_procs=1):
        """
        Execute the binary and return a reference to the subprocess
        created

        :param use_mpi: If mpirun will be called to execute the binary
        :param omp_max_threads: Number of OpenMP threads to be created
                                by default the environment variable
                                OMP_NUM_THREADS is not changed
        :param mpi_num_procs: Number of MPI processes
        :return:
        """
        cwd = os.getcwd()
        os.chdir(self.workdir)
        if omp_max_threads > 0:
            os.environ["OMP_NUM_THREADS"] = str(omp_max_threads)
        if self.stdin_filename is not None:
            self.stdin_file = open(self.stdin_filename, 'r')
        if self.stdout_filename is not None:
            self.stdout_file = open(self.stdout_filename, 'w')
        if self.stderr_filename is not None:
            self.stderr_file = open(self.stderr_filename, 'w')
        if use_mpi:
            sp = subprocess.Popen(['mpirun', '-n', str(mpi_num_procs), self.binary],
                                  stdout=self.stdout_file, stderr=self.stderr_file, stdin=self.stdin_file)
        else:
            sp = subprocess.Popen(self.binary, stdout=self.stdout_file, stderr=self.stderr_file, stdin=self.stdin_file)
        os.chdir(cwd)
        self.runner = sp
        return sp


class CodeInput:

    __metaclass__ = ABCMeta

    def __init__(self):
        self.input_file = None
        self.variables = {}

    @abstractmethod
    def read_inputfile(self):
        pass

    @abstractmethod
    def write_inputfile(self):
        pass

    @abstractmethod
    def __str__(self):
        pass

    def has_variable(self, varname, section=None):
        if self.is_hierarchical:
            if section is None:
                raise ValueError('ERROR: Input variables are hierachical and not section was declared')
            else:
                if section in self.variables and varname in self.variables[section]:
                    return True
                else:
                    return False
        else:
            if varname in self.variables:
                return True
            else:
                return False

    def get_variable(self, varname, section=None):

        if self.is_hierarchical:
            if section is None:
                raise ValueError('ERROR: Input variables are hierachical and not section was declared')
            else:
                if self.has_variable(varname, section=section):
                    return self.variables[section][varname]
                else:
                    return None
        else:
            if self.has_variable(varname):
                return self.variables[varname]
            else:
                return None

    def set_variable(self, varname, value, section=None):
        if self.is_hierarchical:
            if section is None:
                raise ValueError('ERROR: Input variables are hierachical and not section was declared')
            else:
                if self.has_variable(varname, section=section):
                    return self.variables[section][varname]
                else:
                    return None
        else:
            if self.has_variable(varname):
                return self.variables[varname]
            else:
                return None

    @property
    def get_number_variables(self):
        if not self.is_hierarchical:
            return len(self.variables)
        else:
            ret={}
            for i in self.variables:
                ret[i] = len(self.variables[i])
            return ret

    @property
    def is_hierarchical(self):
        return False



class CodeOutput:

    __metaclass__ = ABCMeta

    def __init__(self):
        self.output_file = None
        self.output_values = {}
        self.number_outputs = 0

    @abstractmethod
    def read_outputfile(self):
        pass

    @abstractproperty
    def is_loaded(self):
        pass
