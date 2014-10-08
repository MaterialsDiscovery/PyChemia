__author__ = 'Guillermo Avendano-Franco'

import os as _os
import shutil as _shutil
import subprocess as _subprocess
from multiprocessing import Process
import pychemia


class Runner():

    def __init__(self, code, environment, options):

        if code.lower() not in ['abinit', 'vasp', 'octopus', 'dftbplus', 'fireball']:
            raise ValueError('Code not supported: ', code)
        else:
            self.code = code.lower()

        if environment.lower() not in ['local', 'remote']:
            raise ValueError('Environment must be local or remote: ')
        else:
            self.environment = environment.lower()

        assert (isinstance(options, dict))
        self.options = options
        self.mpi = True
        self.nproc = 2
        if 'mpi' in options:
            self.mpi = options['mpi']
        if 'nproc' in options:
            self.nproc = options['nproc']
        if 'code_bin' in options:
            self.code_bin = options['code_bin']

    def initialize(self, dirpath):
        """
        Utility that copy a given script and execute the given
        command inside the directory
        """
        if not _os.path.isdir(dirpath):
            _os.mkdir(dirpath)
        if self.code == 'abinit':
            pass
            #pychemia.code.abinit.AbiFiles(basedir=dirpath)
        elif self.code == 'vasp':
            pass

    def run(self, dirpath='.'):

        if self.code == 'abinit':
            outfile = 'abinit.stdout'
            errfile = 'abinit.stderr'
            infile = 'abinit.in'
        elif self.code == 'vasp':
            outfile = 'vasp.stdout'
            errfile = 'vasp.stderr'
            infile = None

        if self.environment == 'local':

            def worker(path):
                cwd = _os.getcwd()
                _os.chdir(path)
                outf = open(outfile, 'w')
                errf = open(errfile, 'w')
                if infile is not None:
                    rf = open(infile, 'r')
                else:
                    rf = None
                if self.mpi:
                    print 'Launching '+self.code.upper()+' using MPI with '+str(self.nproc)+' processes'
                    ret = _subprocess.call(['mpirun', '-np', str(self.nproc), self.code_bin],
                                           stdin=rf, stdout=outf, stderr=errf)
                    if ret == 0:
                        print 'Normal completion of '+self.code.upper()
                    else:
                        print 'Finished with error number:', ret
                else:
                    print 'Launching '+self.code.upper()+' using in serial mode'
                    ret = _subprocess.call([self.code_bin], stdin=rf, stdout=outf, stderr=errf)
                    if ret == 0:
                        print 'Normal completion of '+self.code.upper()
                    else:
                        print 'Finished with error number:', ret

                _os.chdir(cwd)
                errf.close()
                outf.close()

            p = Process(target=worker, args=(dirpath, ))
            p.start()
            p.join()

        return p

