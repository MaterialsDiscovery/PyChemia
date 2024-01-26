import datetime
import os
import shutil
import socket
import subprocess
import xml.etree.ElementTree as ElementTree
from pathlib import Path

class SLURM_Runner:
    """
    Class to create and launch submission scripts for Torque/PBS

    """
    def __init__(self, workdir=None, submission_script='batch.slurm', slurm_params=None, input_files=None, jobname=None, command_line=None, use_symlinks=False):
        """
        Class to create and manage job executions on a SLURM resource manager
        """
        self.submission_script = submission_script
        self.jobid = None
        self.join = None
        self.use_symlinks = use_symlinks
        if command_line is None:
            command_line = ""
        else:
            self.command_line = command_line
        if slurm_params is None:
            self.slurm_params = {}
        else:
            self.slurm_params = slurm_params
        if input_files is None:
            self.input_files = []
        else:
            self.input_files = input_files
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


    def write_submission_script(self):
        """
        Write the submisson script based on the contents of slurm_params and code_params

        """
        wf=open(self.workdir + os.sep + self.submission_script, 'w')

        wf.write("#!/bin/bash\n\n")
        for i in self.slurm_params:
            if len(i)==1:
                wf.write("#SBATCH -%s %s\n" % (i, str(self.slurm_params[i])))
            else:
                wf.write("#SBATCH --%s=%s\n" % (i, str(self.slurm_params[i])))

        wf.write("\ncd $SLURM_SUBMIT_DIR\n\n")

        wf.write("%s\n" % self.command_line)

        wf.close()

    def prepare_files(self):

        if not os.path.isdir(self.workdir):
            os.mkdir(self.workdir)
        for i in self.input_files:
            target = os.path.basename(os.path.abspath(i))
            if os.path.lexists(self.workdir + os.sep + target):
                os.remove(self.workdir + os.sep + target)
            if self.use_symlinks:
                os.symlink(Path(i).resolve(),
                    self.workdir + os.sep + target )
            else:
                shutil.copy2(i,self.workdir)
        self.write_submission_script()


    def __repr__(self):
        """
        Evaluatable representation of the object
        """
        return "%s(workdir=%s, submission_script=%s, slurm_params=%s, input_files=%s, jobname=%s, command_line=%s, use_symlinks=%s)" % \
            (self.__class__.__name__,
             "'%s'" % self.workdir if self.workdir is not None else None,
             "'%s'" % self.submission_script if self.submission_script is not None else None,
             "'%s'" % self.slurm_params if self.slurm_params != {} else None,
             "'%s'" % self.input_files if self.input_files != [] else None,
             "'%s'" % self.jobname if self.input_files is not None else None,
             "'%s'" % self.command_line if self.command_line != "" else None,
             "'%s'" % self.use_symlinks)

    def __str__(self):
        """
        String with the contents os the submission script
        """
        ret = "Work Directory = %s\n" % self.workdir
        ret += "Submission Script = %s\n" % self.submission_script
        ret += "SLURM parameters = %s\n" % self.slurm_params
        ret += "Input Files = %s\n" % self.input_files
        ret += "Job Name = %s\n" % self.jobname
        ret += "Command Line = \n\n%s\n" % self.command_line
        return ret

    def submit(self):
        """
        Launches sbatch for the job
        """
        cwd = os.getcwd()
        # Move into the workdir before executing qsub
        os.chdir(self.workdir)
        
        command_line = "sbatch %s" % self.submission_script

        try:
            stdout = subprocess.check_output(command_line, shell=True)
        except subprocess.CalledProcessError as exc:
            print('[ERROR]: SLURM returned error code: %s' % exc.returncode)
            print('         Command line was: %s ' % exc.cmd)
            if len(exc.output) > 0:
                print('         The output from qsub was:\n %s' % exc.output)

        os.chdir(cwd)
        self.jobid = stdout.decode().strip().split()[-1]
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


