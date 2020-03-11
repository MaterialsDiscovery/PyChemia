from abc import ABCMeta, abstractmethod
import os
import subprocess
import collections
import psutil
import shlex
from pychemia import pcm_log

class CodeRun:
    __metaclass__ = ABCMeta

    def __init__(self, executable, workdir='.', use_mpi=False):
        """
        CodeRun is the superclass defining the operations for running codes either directly or asynchronously via a
        submission script on a cluster.

        :param executable:   Name or Path to the executable if can be found on the $PATH or path to the executable
        :param workdir:  Path to a folder where input and output will be located (Default: '.')
        :param use_mpi:  True if code relies on MPI for execution (Default: False)

        """
        self.stdin_file = None
        self.stdout_file = None
        self.stderr_file = None
        self.stdin_filename = None
        self.stdout_filename = None
        self.stderr_filename = None
        self.input_path = None
        self.input = None
        self.executable = executable
        self.workdir = workdir
        self.use_mpi = use_mpi
        self.runner = None

    @abstractmethod
    def set_inputs(self):
        """
        This method must be implemented by child classes, it should write all the input files and prepare environment
        for execution

        :return: None
        """
        pass

    @abstractmethod
    def get_outputs(self):
        """
        This method must be implemented by child classes, it check the existance of output files and prepare the
        reading and parsing of output data.

        :return: None
        """
        pass

    def run(self, num_threads=None, mpi_num_procs=None, nodefile=None, wait=True, verbose=False):
        """
        Run the executable and return a reference to the subprocess
        created. The execution can set a number of threading variables to control their number.
        If the code uses MPI the variable mpi_num_procs control the number of processes used.

        :param num_threads: Control the number of threads in multithreading applications. There are several environment
                            variables used for that purpose. The default is None, meaning that not variable is set,
                            if an integer is given only OMP_NUM_THREADS is set dynamically for the run if a dictionary
                            of variables is provided with keys OMP_NUM_THREADS, OPENBLAS_NUM_THREADS, GOTO_NUM_THREADS
                            and MKL_NUM_THREADS, the corresponding integers associated to those variables are used to
                            set the number of corresponding threads.

        :param mpi_num_procs: Number of MPI processes, if None and the code uses MPI, the default is to set the maximum
                              number of cores available on the system.

        :param nodefile: If a nodefile is provided, and the code uses MPI, the contents are used to select the number of
                         of processes created. Setting nodefile does not override the use of mpi_num_procs but if
                         mpi_num_procs is None it will assume the number of processes equal to the machines
                         declared on nodefile.

        :param wait: If True the process waits until the command is completed, otherwise a Popen instance is returned

        :param verbose: Print extra information before and after the execution

        :return:
        """
        cwd = os.getcwd()
        pcm_log.debug("Current location=%s, workdir=%s" % (cwd, self.workdir))
        os.chdir(self.workdir)

        env_vars = ''
        if num_threads is not None:
            if isinstance(num_threads, int):
                env_vars = 'OMP_NUM_THREADS='+str(num_threads)
        elif isinstance(num_threads, dict):
            for evar in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'GOTO_NUM_THREADS', 'MKL_NUM_THREADS']:
                if evar in num_threads and isinstance(num_threads[evar], int) and num_threads[evar] > 0:
                    env_vars += evar+'='+str(num_threads[evar])+' '

                subprocess.check_output('export %s=%d' % (evar, num_threads), shell=True)
                envvars=subprocess.check_output('echo $%s' % evar, shell=True)
                print(envvars.decode())

        if self.stdin_filename is not None:
            self.stdin_file = open(self.stdin_filename, 'r')
        if self.stdout_filename is not None:
            self.stdout_file = open(self.stdout_filename, 'w')
        if self.stderr_filename is not None:
            self.stderr_file = open(self.stderr_filename, 'w')

        # Checking availability of executable
        if not os.path.exists(self.executable):

            try:
                which_bin = subprocess.check_output('which %s' % self.executable, shell=True)
            except subprocess.CalledProcessError:
                raise ValueError('ERROR: Executable %s could not be found as an absolute, relative or via the $PATH variable' %
                      self.executable)
                
            exec_path = which_bin.decode('utf8').strip()
        else:
            exec_path = self.executable

        pcm_log.debug("Executable: %s " % exec_path)

        if self.use_mpi:

            np = 2
            if mpi_num_procs is None:
                np = psutil.cpu_count(logical=False)
            elif isinstance(mpi_num_procs, int):
                np = mpi_num_procs
            else:
                print("WARNING: Declared variable mpi_num_procs is not integer, defaulting to 2 process")

            command = 'mpirun -n %d %s' % (np, self.executable)
        else:
            command = "%s" % self.executable

        pcm_log.debug("Running: %s" % command)
        process = subprocess.Popen(shlex.split(command), stdout=subprocess.PIPE, stderr=subprocess.PIPE, stdin=self.stdin_file)
        while True:
            output = process.stdout.readline()
            if output == b'' and process.poll() is not None:
                break
            if output != b'' and process.poll() is not None:
                pcm_log.debug("process.poll() is %s" % process.poll())
                pcm_log.debug(output)
            if output and verbose:
                print(output.decode(), end='')
        rc = process.poll()
        self.runner = process
        os.chdir(cwd)
        return process

#       if wait:
#            self.runner.wait()
#            if verbose:
#                print("Program finished with returncode: %d" % self.runner.returncode)
#        os.chdir(cwd)
#        return self.runner


class CodeInput(collections.abc.MutableMapping):

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


class CodeOutput(collections.abc.Mapping):

    __metaclass__ = ABCMeta

    def __init__(self):
        self.output_values = {}

    @abstractmethod
    def read(self):
        pass

    # True means that the run is just complete
    @property
    @abstractmethod
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
