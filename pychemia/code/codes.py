from abc import ABCMeta, abstractmethod, abstractproperty
import os
import subprocess
import collections


class CodeRun:
    __metaclass__ = ABCMeta

    def __init__(self, workdir='.', binary=None):
        self.stdin_file = None
        self.stdout_file = None
        self.stderr_file = None
        self.stdin_filename = None
        self.stdout_filename = None
        self.stderr_filename = None
        self.input_path = None
        self.input = None
        self.binary = binary
        self.workdir = workdir
        self.use_mpi = False

    @abstractmethod
    def set_inputs(self):
        pass

    @abstractmethod
    def get_outputs(self):
        pass

    def execute(self, omp_num_threads=1, mpi_num_procs=1, wait=True):
        """
        Execute the binary and return a reference to the subprocess
        created

        :param wait: If true the process waits until the command is completed
        :param omp_num_threads: Number of OpenMP threads to be created
                                by default the environment variable
                                OMP_NUM_THREADS is not changed
        :param mpi_num_procs: Number of MPI processes
        :return:
        """
        cwd = os.getcwd()
        os.chdir(self.workdir)
        print("Working at: %s" % os.getcwd())
        if omp_num_threads > 0:
            os.environ["OMP_NUM_THREADS"] = str(omp_num_threads)
        if self.stdin_filename is not None:
            self.stdin_file = open(self.stdin_filename, 'r')
        if self.stdout_filename is not None:
            self.stdout_file = open(self.stdout_filename, 'w')
        if self.stderr_filename is not None:
            self.stderr_file = open(self.stderr_filename, 'w')
        if self.use_mpi:
            which_bin = subprocess.check_output('which %s' % self.binary, shell=True)
            print("Executable: %s" % which_bin.decode('utf8').strip())
            command_line = 'mpirun -n %d %s' % (mpi_num_procs, self.binary)
            print("Running: %s" % command_line)
            sp = subprocess.Popen(command_line, shell=True,
                                  stdout=self.stdout_file, 
                                  stderr=self.stderr_file, 
                                  stdin=self.stdin_file)
        else:
            sp = subprocess.Popen("%s" % self.binary, shell=True, 
                                  stdout=self.stdout_file, 
                                  stderr=self.stderr_file, 
                                  stdin=self.stdin_file)
        if wait:
            sp.wait()
            print("Program finished with returncode: %d" % sp.returncode)
        os.chdir(cwd)
        return sp


class CodeInput(collections.MutableMapping):

    __metaclass__ = ABCMeta

    def __init__(self):
        self.input_file = None
        self.variables = {}

    def __contains__(self, x):
        if not self.is_hierarchical:
            return x in self.variables
        else:
            return x[1] in self.variables[x[0]]

    def __delitem__(self, x):
        if not self.is_hierarchical:
            return self.variables.__delitem__(x)
        else:
            return self.variables[x[0]].__delitem__(x[1])

    def __setitem__(self, x, value):
        if not self.is_hierarchical:
            return self.variables.__setitem__(x, value)
        else:
            return self.variables[x[0]].__setitem__(x[1], value)

    def __getitem__(self, x):
        if not self.is_hierarchical:
            return self.variables.__getitem__(x)
        else:
            return self.variables[x[0]].__getitem__(x[1])

    def __iter__(self):
        return self.variables.__iter__()

    def __len__(self):
        return self.variables.__len__()

    @abstractmethod
    def read(self):
        pass

    @abstractmethod
    def __str__(self):
        pass

    def write(self, filename=None):
        """
        Write an input object into a text
        file that ABINIT can use as an input
        file

        Args:
            filename:
                The 'abinit.in' filename that will be written
        """
        if filename is None:
            if self.input_file is not None:
                filename = self.input_file
            else:
                raise ValueError("Not filename indicated")
        wf = open(filename, 'w')
        wf.write(self.__str__())
        wf.close()

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
                raise ValueError('ERROR: Input variables are hierarchical and not section was declared')
            else:
                self.variables[section][varname] = value
        else:
            self.variables[varname] = value

    @property
    def get_number_variables(self):
        if not self.is_hierarchical:
            return len(self.variables)
        else:
            ret = {}
            for i in self.variables:
                ret[i] = len(self.variables[i])
            return ret

    @property
    def is_hierarchical(self):
        return False


class CodeOutput(collections.Mapping):

    __metaclass__ = ABCMeta

    def __init__(self):
        self.output_values = {}

    @abstractmethod
    def read(self):
        pass

    # True means that the run is just complete
    @abstractproperty
    def is_finished(self):
        pass

    @property
    def is_loaded(self):
        return not self.output_values == {}    

    def __contains__(self, x):
        return x in self.output_values

    def __getitem__(self, x):
        return self.output_values.__getitem__(x)
    
    def __iter__(self):
        return self.output_values.__iter__()

    def __len__(self):
        return self.output_values.__len__()
