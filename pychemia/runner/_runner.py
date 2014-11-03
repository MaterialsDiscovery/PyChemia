__author__ = 'Guillermo Avendano-Franco'

import os
import subprocess as _subprocess
from multiprocessing import Process
import time


class Runner():

    def __init__(self, code, code_bin, environment, use_mpi=True, nmpiproc=2, nconcurrent=1):

        if code.lower() not in ['abinit', 'vasp', 'octopus', 'dftbplus', 'fireball']:
            raise ValueError('Code not supported: ', code)
        else:
            self.code = code.lower()

        if environment.lower() not in ['local', 'remote']:
            raise ValueError('Environment must be local or remote: ')
        else:
            self.environment = environment.lower()

        self.use_mpi = use_mpi
        self.nmpiproc = nmpiproc
        self.nconcurrent = nconcurrent
        self.code_bin = code_bin

    def initialize(self, dirpath):
        """
        Utility that copy a given script and execute the given
        command inside the directory
        """
        if not os.path.isdir(dirpath):
            os.mkdir(dirpath)
        if self.code == 'abinit':
            pass
            #pychemia.code.abinit.AbiFiles(basedir=dirpath)
        elif self.code == 'vasp':
            pass

    def run(self, dirpath='.', verbose=False):

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
                cwd = os.getcwd()
                os.chdir(path)
                outf = open(outfile, 'w')
                errf = open(errfile, 'w')
                if infile is not None:
                    rf = open(infile, 'r')
                else:
                    rf = None
                if self.use_mpi:
                    #print 'Launching '+self.code.upper()+' using MPI with '+str(self.nmpiproc)+' processes'
                    ret = _subprocess.call(['mpirun', '-np', str(self.nmpiproc), self.code_bin],
                                           stdin=rf, stdout=outf, stderr=errf)
                else:
                    #print 'Launching '+self.code.upper()+' using in serial mode'
                    ret = _subprocess.call([self.code_bin], stdin=rf, stdout=outf, stderr=errf)

                if verbose:
                    if ret == 0:
                        print 'Normal completion of '+self.code.upper()
                    else:
                        print 'Finished with error number:', ret

                os.chdir(cwd)
                errf.close()
                outf.close()

            p = Process(target=worker, args=(dirpath, ))
            p.start()
            p.join()

        return p

    def run_multidirs(self, workdirs, worker, checker):
        pt = []
        for i in range(self.nconcurrent):
            pt.append(None)

        index = 0
        while True:
            for i in range(self.nconcurrent):
                if pt[i] is None or not pt[i].is_alive():
                    ret = checker(workdirs[index])
                    if ret:
                        print 'Launching for process '+str(i)+' on dir '+workdirs[index]+' index '+str(index)
                        pt[i] = Process(target=worker, args=(workdirs[index],))
                        pt[i].start()
                    index = (index+1) % len(workdirs)
            time.sleep(30)
            complete = True
            for idir in workdirs:
                if not os.path.isfile(idir+os.sep+'COMPLETE'):
                    complete = False
            if complete:
                print 'Finishing...'
                break
