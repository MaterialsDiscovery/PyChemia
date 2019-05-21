#!/usr/bin/env python

import os                              # Operating System module (create directories)
import time                            # To create a sleeping timer
import random
import subprocess                      # To launch an external command (Such as maise)
from multiprocessing import Process    # Multiprocessing python (Concurrency)


def maise_worker(case, worker_args):
    # os.chdir(workdir)
    # subprocess.call(['maise'])
    random.seed()
    workdir = "MAISE_%03d" % case
    if not os.path.isdir(workdir):
        os.mkdir(workdir)
    wf = open(workdir+'/EXECTUING MAISE', 'w')
    wf.write('Put results here')
    ttl = random.randint(10, 40)
    print('Evaluating case: %d this calculation will take %d seconds' % (case, ttl))
    time.sleep(ttl)
    wf.close()


def run(to_evaluate, worker, nconcurrent=1, sleeping_time=120):
    """
    For a given list of directories to evaluate, start processing concurrently the worker until all directories are
    processed

    :param worker:        Function to evaluate for each case
    :param sleeping_time: Number of seconds to wait until checking for complete jobs
    :param to_evaluate: (list) List of directories to evaluate
    :param nconcurrent: (int)  number of concurrent evaluations of the worker
    """

    procs = []
    worker_args = (None,)  # This is a dummy argument, not playing any role here, it must be a tuple always

    # Creating a list to store the 'nconcurrent' running jobs
    for i in range(nconcurrent):
        procs.append(None)

    # Main loop looking permanently for candidates for evaluation
    while True:

        index = 0
        currently_evaluating = 0
        for j in range(nconcurrent):
            if procs[j] is not None and procs[j].is_alive():
                currently_evaluating += 1
        print('Candidates to evaluate: %d  Candidates in evaluation: %d' % (len(to_evaluate), currently_evaluating))

        while index < len(to_evaluate):

            slot = None
            while True:
                for j in range(nconcurrent):
                    if procs[j] is None or not procs[j].is_alive():
                        slot = j
                        break
                if slot is None:
                    time.sleep(sleeping_time)
                else:
                    break

            case = to_evaluate[index]
            print("I will run case: %d with index: %d on slot %d/%d" % (case, index, slot+1, nconcurrent))

            # This is the actual call to the worker, it must be a function with 4 arguments:
            # The database settings, the entry identifier, the working directory and arguments for the worker
            procs[slot] = Process(target=worker, args=(case, worker_args))
            procs[slot].start()
            time.sleep(1)

            index += 1
        time.sleep(sleeping_time)

        # Condition to finish execution, if removed the code will run forever looking for more cases to run
        if index >= len(to_evaluate):
            break


if __name__ == '__main__':

    to_evaluate = range(100)  # Creating 100 cases to evaluate
    nconcurrent = 16          # Number of concurrent executions

    run(to_evaluate=to_evaluate, nconcurrent=nconcurrent, worker=maise_worker, sleeping_time=5)


