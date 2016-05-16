import os
import subprocess
from multiprocessing import Process
import time


class Runner:
    def __init__(self, code, code_bin, environment, use_mpi=True, nmpiproc=2, nconcurrent=1, runtime=3600):

        if code.lower() not in ['abinit', 'vasp', 'octopus', 'dftb', 'fireball']:
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
        self.runtime = runtime

    def initialize(self, dirpath):
        """
        Utility that copy a given script and execute the given
        command inside the directory

        :param dirpath: (str) Directory to execute runner
        """
        if not os.path.isdir(dirpath):
            os.mkdir(dirpath)
        if self.code == 'abinit':
            pass
            # pychemia.code.abinit.AbiFiles(basedir=dirpath)
        elif self.code == 'vasp':
            pass

    def run(self, dirpath='.', analyser=None):

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
                outf = open(outfile, 'a')
                errf = open(errfile, 'a')
                for i in [outf, errf]:
                    i.write('' + (40 * '=') + ' New Run ' + (40 * '=') + '\n')
                if infile is not None:
                    rf = open(infile, 'r')
                else:
                    rf = None

                initime = time.time()
                if self.use_mpi:
                    childp = subprocess.Popen(['mpirun', '-np', str(self.nmpiproc),
                                               '--map-by', 'socket:PE=2', self.code_bin],
                                              stdin=rf, stdout=outf, stderr=errf)
                else:
                    childp = subprocess.Popen([self.code_bin], stdin=rf, stdout=outf, stderr=errf)

                childp.poll()
                while childp.returncode is None:
                    childp.poll()
                    rf = open('vasp.stdout', 'r')
                    for iline in rf.readlines()[-20:]:
                        if iline.strip().startswith("WARNING: Sub-Space-Matrix is not hermitian in DAV"):
                            print("[WARNING: Sub-Space-Matrix is not hermitian in DAV] Stopping now...")
                            self._stop_run(childp)
                    if time.time() - initime > self.runtime:
                        self._stop_run(childp)
                        break
                    else:
                        if analyser is not None:
                            ret = analyser()
                            print(os.path.basename(path), ret)
                        time.sleep(60)

                print('The return code was', childp.returncode)
                os.chdir(cwd)
                errf.close()
                outf.close()

            p = Process(target=worker, args=(dirpath,))
            p.start()
            p.join()

        return p

    @staticmethod
    def _stop_run(childp, softtime=600, extratime=60):
        print('Creating Stoping files')
        wf = open('CHKPT', 'w')
        wf.close()
        wf = open('STOPCAR', 'w')
        wf.write('LSTOP = .TRUE.\n')
        wf.close()
        time.sleep(softtime)

        childp.poll()
        if childp.returncode is None:
            print('Sending SIGTERM to process')
            childp.terminate()
            time.sleep(extratime)

        childp.poll()
        if childp.returncode is None:
            print('Sending SIGKILL to process')
            childp.kill()
        if os.path.isfile('CHKPT'):
            os.remove('CHKPT')
        if os.path.isfile('STOPCAR'):
            os.remove('STOPCAR')

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
                        print('Launching for process ' + str(i) + ' on dir ' + os.path.basename(
                            workdirs[index]) + ' index ' + str(index))
                        pt[i] = Process(target=worker, args=(workdirs[index],))
                        pt[i].start()
                    index = (index + 1) % len(workdirs)
            time.sleep(30)
            complete = True
            for idir in workdirs:
                if not os.path.isfile(idir + os.sep + 'COMPLETE'):
                    complete = False
            if complete:
                print('Finishing...')
                break

    def run_multidirs_nonstop(self, workdirs, worker, checker):

        def superworker(s_workdirs, s_worker, s_checker):
            self.run_multidirs(s_workdirs, s_worker, s_checker)

        proc = Process(target=superworker, args=(workdirs, worker, checker))
        proc.start()
