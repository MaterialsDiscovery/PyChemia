import datetime
import os
import socket
import subprocess
import xml.etree.ElementTree as ElementTree


class PBSRunner:
    def __init__(self, workdir, filename='batch.pbs'):

        self.template = None
        self.walltime = None
        self.template = None
        self.queue = None
        self.mail = None
        self.message = None
        self.walltime = None
        self.ppn = None
        self.nodes = None
        self.workdir = workdir
        self.filename = filename
        if self.workdir[-1] == os.sep:
            self.workdir = self.workdir[:-1]
        self.name = os.path.basename(self.workdir)

    def initialize(self, nodes=1, ppn=2, walltime=None, message='ae', mail=None, queue=None):

        if walltime is None:
            walltime = [12, 0, 0]
        self.set_walltime(walltime)
        self.nodes = nodes
        self.ppn = ppn
        self.message = message
        self.mail = mail
        self.queue = queue

    def set_walltime(self, walltime):

        if len(walltime) == 1:
            walltime = [0] + walltime
        if len(walltime) == 2:
            walltime = [0] + walltime
        if len(walltime) == 3:
            walltime = [0] + walltime

        self.walltime = walltime

    def set_template(self, template):
        if os.path.isfile(template):
            self.template = open(template).read()
        elif isinstance(template, str):
            self.template = template

    def write_pbs(self):

        wf = open(self.workdir + os.sep + self.filename, 'w')
        wt = self.walltime
        wf.write("""#!/bin/sh

#PBS -N %s
#PBS -l nodes=%d:ppn=%d
#PBS -l walltime=%d:%02d:%02d
#PBS -m %s
#PBS -k n
""" % (self.name, self.nodes, self.ppn, wt[0] * 24 + wt[1], wt[2], wt[3], self.message))

        if self.mail is not None:
            wf.write("#PBS -M %s\n" % self.mail)
        if self.queue is not None:
            wf.write("#PBS -q %s\n" % self.queue)
        wf.write('\ncd $PBS_O_WORKDIR\n')

        if self.template is not None:
            wf.write("%s\n" % self.template)
        wf.close()

    def submit(self, priority=0):
        cwd = os.getcwd()
        os.chdir(self.workdir)
        returncode = subprocess.call(["qsub", "%s" % self.filename, '-p', '%d' % priority])
        if returncode != 0:
            print('Some error happended:', returncode)
        os.chdir(cwd)


def get_jobs(user):
    data = subprocess.check_output(['qstat', '-x', '-f', '-u', user])
    xmldata = ElementTree.fromstring(data)
    jobs = [i.find('Job_Name').text for i in xmldata.findall('Job')]
    return jobs


def report_cover():
    ret = 'PyChemia Execution Report\n'
    ret += '=========================\n\n'
    ret += 'Hostname\n'
    ret += '--------\n\n'
    ret += socket.gethostname() + '\n\n'
    ret += 'Date\n'
    ret += '----\n'
    dt = datetime.datetime.now()
    ret += dt.strftime("%A, %d. %B %Y %I:%M%p") + '\n\n'

    ret += 'PBS VARIABLES\n'
    ret += '-------------\n\n'
    for x in ['PBS_O_HOST', 'PBS_SERVER', 'PBS_O_QUEUE', 'PBS_O_WORKDIR', 'PBS_ARRAYID', 'PBS_ENVIRONMENT',
              'PBS_JOBID', 'PBS_JOBNAME', 'PBS_NODEFILE', 'PBS_QUEUE']:
        if os.getenv(x) is not None:
            ret += x + ' = ' + os.getenv(x) + '\n'
    return ret
