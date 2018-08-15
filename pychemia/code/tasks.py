import os
import json
import subprocess
from abc import ABCMeta, abstractmethod


class Task:
    __metaclass__ = ABCMeta

    def __init__(self, structure, task_params=None, workdir='.', executable=None):
        self.structure = structure
        self.workdir = workdir
        if not os.path.isdir(self.workdir):
            os.mkdir(self.workdir)

        if task_params is not None:
            self.task_params = task_params
        else:
            self.task_params = {}
        self.executable = executable
        self.finished = False
        self.success = False
        self.started = False
        self.output = {}
        self.report_dir = self.workdir + os.sep + 'REPORT'

    @abstractmethod
    def run(self):
        pass

    @abstractmethod
    def plot(self):
        pass

    @abstractmethod
    def report(self):
        pass

    @abstractmethod
    def load(self):
        pass

    def save(self, filename=None):
        if filename is None:
            filename = self.workdir + os.sep + 'task.json'
        wf = open(filename, 'w')
        ret = {'task_params': self.task_params, 'output': self.output}
        json.dump(ret, wf, sort_keys=True, indent=4, separators=(',', ': '))
        wf.close()

    def status(self):
        if self.finished:
            print('Task finished')
        if self.started:
            print('Task started')
        if self.success:
            print('Task completed successfully')

    def report_end(self, html, file_format):
        from lxml import etree

        doctype = '<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.1//EN" "http://www.w3.org/TR/xhtml11/DTD/xhtml11.dtd">'
        result = etree.tostring(html,
                                xml_declaration=True,
                                doctype=doctype,
                                encoding='utf-8',
                                standalone=False,
                                with_tail=False,
                                method='xml',
                                pretty_print=True)

        wf = open(self.report_dir + os.sep + 'index.html', 'w')
        wf.write(result)

        if file_format != 'html':
            cwd = os.getcwd()
            os.chdir(self.report_dir)
            stderr = open('pandoc_out.log', 'w')
            stdout = open('pandoc_err.log', 'w')
            sp = subprocess.Popen(['pandoc', 'index.html', '-o', 'visual.' + file_format], stderr=stderr, stdout=stdout)
            os.chdir(cwd)
            stderr.close()
            stdout.close()
            return sp
