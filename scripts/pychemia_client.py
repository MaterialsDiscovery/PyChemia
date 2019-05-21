#!/usr/bin/env python

import getopt
import os
import random
import signal
import socket
import subprocess
import sys
import time

from pychemia.utils.computing import get_int

__author__ = 'Guillermo Avendano-Franco'


class TimedOutExc(Exception):
    pass


def deadline(timeout, *args):
    def decorate(f):
        def handler(signum, frame):
            print('signum:', signum)
            print('frame:', frame)
            raise TimedOutExc()

        def new_f(*args):
            signal.signal(signal.SIGALRM, handler)
            signal.alarm(timeout)
            return f(*args)

        new_f.__name__ = f.__name__
        return new_f

    return decorate


def usage(name):
    print("""
NAME
    %s

DESCRIPTION
    PyChemia Client intended to be run by compute nodes


OPTIONS

    --help, -h
        Return information on the options and use of this script

    --port, -p <int> (Default: 10000)
        The UDP port that will be used for transfer data with clients

    --workdir, w <string>
        Working directory were files will be created and jobs runned

    --ip, -i <string>
        IP address to bind the server

""" % os.path.basename(name))


@deadline(60)
def inquirer(ip, port):
    # Create a UDP socket
    sock = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
    server_address = (ip, port)
    print('Connecting to %s:%d' % (ip, port), file=sys.stderr)

    try:
        message = 'COUNT'
        # Send data
        print('sending "%s"' % message, file=sys.stderr)
        sock.sendto(message, server_address)

        # Receive response
        print('waiting to receive', file=sys.stderr)
        data, server = sock.recvfrom(4096)
        print('received "%s"' % data, file=sys.stderr)

        message = 'WORKDIR'
        # Send data
        print('sending "%s"' % message, file=sys.stderr)
        sock.sendto(message, server_address)

        # Receive response
        print('waiting to receive', file=sys.stderr)
        data, server = sock.recvfrom(4096)
        print('received "%s"' % data, file=sys.stderr)
        workdir = data

        message = 'GET'
        # Send data
        print('sending "%s"' % message, file=sys.stderr)
        sock.sendto(message, server_address)

        # Receive response
        print('waiting to receive', file=sys.stderr)
        data, server = sock.recvfrom(4096)
        print('received "%s"' % data, file=sys.stderr)
        entry_id = data

    finally:
        print('closing socket', file=sys.stderr)
        sock.close()

    return workdir + os.sep + entry_id, entry_id


def executor(workdir):
    print('Work directory', workdir)
    os.chdir(workdir)
    subprocess.call(['python', 'executor.py'])


@deadline(60)
def finisher(entry_id, ip, port):
    # Create a UDP socket
    sock = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)

    server_address = (ip, port)

    try:
        message = 'FINISHED:' + str(entry_id)
        # Send data
        print('sending "%s"' % message, file=sys.stderr)
        sock.sendto(message, server_address)

        # Receive response
        print('waiting to receive', file=sys.stderr)
        data, server = sock.recvfrom(4096)
        print('received "%s"' % data, file=sys.stderr)

    finally:
        print('closing socket', file=sys.stderr)
        sock.close()


def main(argv):
    try:
        opts, args = getopt.getopt(argv[1:], "hi:p:", ["help", "ip=", "port="])
    except getopt.GetoptError:
        print('Wrong parsing of options')
        usage(argv[0])
        sys.exit(2)

    # Default Values
    port = None
    ip = None

    for opt, arg in opts:

        if opt in ("-h", "--help"):
            usage(argv[0])
            sys.exit()
        elif opt in ("-p", "--port"):
            port = get_int(arg)
        elif opt in ("-i", "--ip"):
            ip = arg

    if ip is None:
        ip = socket.gethostbyname(socket.gethostname())

    print('PyChemia Client connecting to server on : ', ip)

    workdir = None
    entry_id = None
    while True:
        if port is not None:
            try:
                workdir, entry_id = inquirer(ip, port)
            except TimedOutExc as e:
                print("INQUIRER took too long: %s", e.message)
        else:
            lt = time.localtime()
            random.seed(lt.tm_yday * 24 + lt.tm_hour)
            rndport = random.randint(10000, 20000)
            print('PyChemia Client using port : ', rndport)
            try:
                workdir, entry_id = inquirer(ip, rndport)
            except TimedOutExc as e:
                print("Error({0}): {1}".format(e, e.message))

        print('We got the workdir: ', workdir)
        print('The entry_id is : ', entry_id)

        if workdir is not None and entry_id is not None:
            if entry_id == '':
                print('No entries to execute, taking a nap')
                time.sleep(15)
            else:
                break
        else:
            'Inquierer failed to get data, retry...'

    executor(workdir)

    if port is not None:
        finisher(entry_id, ip, port)
    else:
        lt = time.localtime()
        random.seed(lt.tm_yday * 24 + lt.tm_hour)
        rndport = random.randint(10000, 20000)
        print('PyChemia Client using port : ', rndport)

        try:
            finisher(entry_id, ip, rndport)
        except TimedOutExc as e:
            print("Error({0}): {1}".format(e, e.message))


if __name__ == "__main__":
    main(sys.argv)
