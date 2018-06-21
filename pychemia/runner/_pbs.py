import datetime
import os
import socket
import subprocess
import xml.etree.ElementTree as ElementTree


class PBSRunner:
    """
    Class to create and launch submission scripts for Torque/PBS

    """
    def __init__(self, workdir=None, filename='batch.pbs', jobname=None, template=None):
        """
        Class to create and manage job executions on a PBS queue system
        """
        self.template = None
        self.walltime = None
        self.template = None
        self.template_text = None
        self.queue = None
        self.mail = None
        self.message = None
        self.walltime = None
        self.ppn = None
        self.nodes = None
        self.features = None
        self.filename = filename
        self.jobid = None
        self.pvmem = None
        self.join = None
        if workdir is None:
            self.workdir = "."
        elif workdir[-1] == os.sep:
            self.workdir = workdir[:-1]
        else:
            self.workdir = workdir
        if jobname is None:
            self.jobname = os.path.basename(os.path.abspath(self.workdir))
        else:
            self.jobname = jobname
        if template is not None:
            self.set_template(template)
            self.template = template

    def set_pbs_params(self, nodes=None, ppn=None, walltime=None, message=None, mail=None, queue=None, features=None,
                       join=None, pvmem=None):
        """
        Set various PBS params

        """
        if walltime is None:
            walltime = [12, 0, 0]
        self.set_walltime(walltime)
        if nodes is not None: 
            self.nodes = nodes
        if self.nodes is None:
            self.nodes = 1
        self.ppn = ppn
        self.message = message
        self.mail = mail
        self.queue = queue
        self.features = features
        self.pvmem = pvmem
        if join is not None:
            if join in ['eo', 'oe']:
                self.join = join
            else:
                self.join = None
        else:
            self.join = None

    def set_walltime(self, walltime):
        """
        Ensures a list of 4 entries for the walltime
        the variable must be [days, hours, minutes, seconds]
        if the list contains less than 4 values completes the others with zeros
        """

        if len(walltime) == 1:
            walltime = [0] + walltime
        if len(walltime) == 2:
            walltime = [0] + walltime
        if len(walltime) == 3:
            walltime = [0] + walltime
        self.walltime = walltime

    def set_template(self, template):
        """
        Prepares the template text, the actual code to be executed on the cluster.

        """
        if os.path.isfile(template):
            self.template_text = open(template).read()
        elif isinstance(template, str):
            self.template_text = template

    def __repr__(self):
        """
        Evaluatable representation of the object
        """
        return "%s(workdir=%s, filename=%s, jobname=%s, template=%s)" % \
            (self.__class__.__name__,
             "'%s'" % self.workdir if self.workdir is not None else None,
             "'%s'" % self.filename if self.filename is not None else None,
             "'%s'" % self.jobname if self.jobname is not None else None,
             "'%s'" % self.template if self.template is not None else None)

    def __str__(self):
        """
        String with the contents os the submission script
        """
        wt = self.walltime
        if self.features is None:
            feat = ':'
        else:
            feat = ':'+self.features+':'
        if self.nodes is None:
            self.nodes = 1

        ret = "#!/bin/sh\n"
        if self.jobname is not None:
            ret += "\n#PBS -N %s\n" % self.jobname
        if self.nodes is not None:
            ret += "#PBS -l nodes=%d" % self.nodes
        if self.ppn is not None:
            ret += "%sppn=%d\n" % (feat, self.ppn)
        if self.walltime is not None:
            ret += "#PBS -l walltime=%d:%02d:%02d\n" % (wt[0] * 24 + wt[1], wt[2], wt[3])
        if self.pvmem is not None:
            ret += "#PBS -l pvmem=%s\n" % self.pvmem 
        if self.message is not None:
            ret += "#PBS -m %s\n" % self.message
        if self.mail is not None:
            ret += "#PBS -M %s\n" % self.mail
        if self.queue is not None:
            ret += "#PBS -q %s\n" % self.queue
        if self.join is not None:
            ret += "#PBS -j %s\n" % self.join 

        ret += '\ncd $PBS_O_WORKDIR\n'

        if self.template is not None:
            ret += "\n# from template\n"
            ret += "%s\n" % self.template_text
        return ret

    def write(self):
        """
        Writes the submission script on file
        """
        if not os.path.isdir(self.workdir):
            os.mkdir(self.workdir)
        wf = open(self.workdir + os.sep + self.filename, 'w')
        wf.write(str(self))
        wf.close()

    def submit(self, priority=None):
        """
        Launches qsub for the job
        """
        cwd = os.getcwd()
        if not os.path.isdir(self.workdir):
            os.mkdir(self.workdir)
        if not os.path.isfile(self.workdir + os.sep + self.filename):
            self.write()

        # Move into the workdir before executing qsub
        os.chdir(self.workdir)
        if priority is None:
            command_line = "qsub %s" % self.filename
        else:
            command_line = "qsub %s -p %d" % (self.filename, priority)

        try:
            stdout = subprocess.check_output(command_line, shell=True)
        except subprocess.CalledProcessError as exc:
            print('[ERROR]: Torque/PBS returned error code: %s' % exc.returncode)
            print('         Command line was: %s ' % exc.cmd)
            if len(exc.output) > 0:
                print('         The output from qsub was:\n %s' % exc.output)

        os.chdir(cwd)
        self.jobid = stdout.decode().strip()
        return stdout.decode().strip()

    def get_stats(self, username):
        if self.jobid is None:
            return None
        else:
            jobs = get_jobs(username)
            if self.jobid in jobs:
                return jobs[self.jobid]
        

def get_jobs(user):
    """
    Check all the jobs submitted by a given 'user'
    Returns a dictionary with the JobIDs and names.
    """
    data = subprocess.check_output(['qstat', '-x', '-f', '-u', user])
    xmldata = ElementTree.fromstring(data)
    jobs = xmldata.findall('Job')
    ret = {}
    for ijob in jobs:
        children = ijob.getchildren()
        jobid = ijob.findall('Job_Id')[0].text
        ret[jobid] = {}
        for child in children:
            ret[jobid][child.tag] = child.text
    return ret


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
