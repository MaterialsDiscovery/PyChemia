__author__ = 'macbook'

import os
import subprocess


class PBSRunner():

    def __init__(self, workdir, filename='batch.pbs'):

        self.template = template
        self.walltime = walltime
        self.template = None
        self.queue = queue
        self.mail = mail
        self.message = message
        self.walltime = None
        self.ppn = ppn
        self.nodes = nodes
        self.name = os.path.basename(self.workdir)
        self.workdir = workdir
        self.filename = filename
        if self.workdir[-1] == os.sep:
            self.workdir = self.workdir[:-1]

    def initialize(self, name=None, nodes=1, ppn=2, walltime=None, message='ae', mail=None, queue=None):

        if name is None:
            walltime=[12, 0, 0]
        self.set_walltime(walltime)

    def set_walltime(self, walltime):

        if len(walltime) == 1:
            self.walltime = [0] + walltime
        if len(walltime) == 2:
            self.walltime = [0] + walltime
        if len(walltime) == 3:
            self.walltime = [0] + walltime

    def set_template(self, template):
        pass

    def write_pbs(self):

        wf = open(self.workdir+os.sep+self.filename, 'w')
        wt = self.walltime
        wf.write("""#!/bin/sh

#PBS -N %s
#PBS -l nodes=%d:ppn=%d
#PBS -l walltime=%d:%02d:%02d
#PBS -m %s
""" % (self.name, self.nodes, self.ppn, wt[0]*24+wt[1], wt[2], wt[3], self.message))

        if self.mail is not None:
            wf.write("#PBS -M %s\n" % self.mail)
        if self.queue is not None:
            wf.write("#PBS -q %s\n" % self.queue)
        wf.write('\ncd %s\n' % os.path.abspath(self.workdir))

        if self.template is not None:
            wf.write("\n%s\n" % self.template)
        wf.close()

    def submit(self):
        cwd = os.getcwd()
        os.chdir(self.workdir)
        returncode = subprocess.call(["qsub", "%s" % self.filename])
        if returncode != 0:
            print 'Some error happended:', returncode
        os.chdir(cwd)
