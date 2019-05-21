#!/usr/bin/env python

import os
import getopt
import json
import multiprocessing
import random
import socket
import sys
import time
from bson.objectid import ObjectId
import pychemia
from pychemia.utils.computing import get_int

__author__ = 'Guillermo Avendano-Franco'


def usage(name):
    print("""
NAME
    %s

DESCRIPTION
    PyChemia Server intended to be run as a daemon on the frontend


OPTIONS

    --help, -h
        Return information on the options and use of this script

    --dbsettings_file, -d <string>
        The PyChemiaQueue database that will be used as source for jobs

    --port, -p <int> (Default: 10000)
        The UDP port that will be used for transfer data with clients

    --workdir, w <string>
        Working directory were files will be created and jobs runned

    --ip, -i <string>
        IP address to bind the server

""" % os.path.basename(name))


def get_database(dbsettings):
    if 'host' not in dbsettings:
        dbsettings['host'] = 'localhost'
    if 'port' not in dbsettings:
        dbsettings['port'] = 27017
    if 'ssl' not in dbsettings:
        dbsettings['ssl'] = False

    if 'user' not in dbsettings:
        pcq = pychemia.db.PyChemiaQueue(name=dbsettings['name'], host=dbsettings['host'], port=dbsettings['port'],
                                        ssl=dbsettings['ssl'])
    else:
        pcq = pychemia.db.PyChemiaQueue(name=dbsettings['name'], host=dbsettings['host'], port=dbsettings['port'],
                                        user=dbsettings['user'], passwd=dbsettings['passwd'], ssl=dbsettings['ssl'])
    return pcq


def listener(dbsettings, ip, port, workdir):
    # Create a TCP/IP socket
    sock = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)

    # Bind the socket to the port
    server_address = (ip, port)
    print('starting up on %s port %s' % server_address, file=sys.stderr)
    sock.bind(server_address)

    pcq = get_database(dbsettings)

    while True:
        print('\nwaiting to receive message', file=sys.stderr)
        data, address = sock.recvfrom(4096)

        print('received %s bytes from %s' % (len(data), address), file=sys.stderr)
        print(data, file=sys.stderr)

        if data == 'COUNT':
            ans = str(pcq.db.pychemia_entries.count())
            sent = sock.sendto(ans, address)
            print('sent %s bytes back to %s' % (sent, address), file=sys.stderr)
        if data == 'WORKDIR':
            sent = sock.sendto(workdir, address)
            print('sent %s bytes back to %s' % (sent, address), file=sys.stderr)
        if data == 'GET':
            # Selectiong an entry not submitted for execution
            entry = pcq.db.pychemia_entries.find_one({'meta.submitted': False}, {'_id': 1})
            if entry is not None:
                entry_id = entry['_id']
                pcq.db.pychemia_entries.update({'_id': entry_id}, {'$set': {'meta.submitted': True}})
                deploy(entry_id, pcq, workdir)
            else:
                entry_id = ''
            sent = sock.sendto(str(entry_id), address)
            print('sent %s bytes back to %s' % (sent, address), file=sys.stderr)
        if data.startswith('FINISHED:'):
            entry_id = ObjectId(data.split(':')[1])
            collect(entry_id, pcq, workdir)
            pcq.db.pychemia_entries.update({'_id': entry_id}, {'$set': {'meta.finished': True}})
            sent = sock.sendto('OK', address)
            print('sent %s bytes back to %s' % (sent, address), file=sys.stderr)


def deploy(entry_id, pychemia_queue, basedir):
    workdir = os.path.abspath(basedir) + os.sep + str(entry_id)
    if not os.path.exists(workdir):
        os.mkdir(workdir)
    pychemia_queue.write_input_files(entry_id, destination=workdir)

    entry = pychemia_queue.db.pychemia_entries.find_one({'_id': entry_id}, {'job': 1})
    job_settings = entry['job']
    print(job_settings)

    if len(job_settings) > 0:
        wf = open(workdir + os.sep + 'job.json', 'w')
        json.dump(job_settings, wf)
        wf.close()

    st = pychemia_queue.get_input_structure(entry_id)
    st.save_json(workdir + os.sep + 'structure.json')

    inp = pychemia_queue.get_input_variables(entry_id)
    if inp is not None:
        wf = open(workdir + os.sep + 'input.json', 'w')
        json.dump(inp, wf)
        wf.close()


def collect(entry_id, pychemia_queue, workdir):
    destination = os.path.abspath(workdir) + os.sep + str(entry_id)
    if os.path.isfile(destination + os.sep + 'results.json'):
        results = json.load(open(destination + os.sep + 'results.json'))
        pychemia_queue.db.pychemia_entries.update({'_id': entry_id}, {'$set': {'output': results}})


def main(argv):
    try:
        opts, args = getopt.getopt(argv[1:], "hd:p:i:w:", ["help", "dbsettings_file=", "port=", "ip=", 'workdir='])
    except getopt.GetoptError:
        print('Error with options')
        usage(argv[0])
        sys.exit(2)

    if len(opts) == 0:
        usage(argv[0])
        sys.exit(2)

    # Default Values
    workdir = None
    dbsettings_file = None
    port = None
    ip = None

    for opt, arg in opts:

        if opt in ("-h", "--help"):
            usage(argv[0])
            sys.exit()
        elif opt in ("-w", "--workdir"):
            workdir = arg
        elif opt in ("-d", "--dbsettings_file"):
            dbsettings_file = arg
        elif opt in ("-p", "--port"):
            port = get_int(arg)
        elif opt in ("-i", "--ip"):
            ip = arg

    if workdir is None:
        print('Missing workdir')
        usage(argv[0])
        sys.exit()
    if not os.path.isdir(workdir):
        os.mkdir(workdir)

    print('dbsettings_file:', dbsettings_file)
    if dbsettings_file is None or not os.path.exists(dbsettings_file):
        print('Missing dbsettings')
        usage(argv[0])
        sys.exit()
    else:
        dbsettings = json.load(open(dbsettings_file))

    if ip is None:
        ip = socket.gethostbyname(socket.gethostname())

    print('PyChemia Server binded to IP: ', ip)

    pcq = get_database(dbsettings)
    print('Connected with database, number of entries: ', pcq.db.pychemia_entries.count())

    if port is not None:
        print('Using port : %s' % port)
        p = multiprocessing.Process(target=listener, args=(dbsettings, ip, port, workdir))
        p.start()
    else:
        old_port = 0
        p = None
        while True:
            lt = time.localtime()
            random.seed(lt.tm_yday * 24 + lt.tm_hour)
            port = random.randint(10000, 20000)
            if port != old_port:
                if p is not None and p.is_alive():
                    print('Terminating listener on port : %s' % old_port)
                    p.terminate()
                print('Using port : %s' % port)
                p = multiprocessing.Process(target=listener, args=(dbsettings, ip, port, workdir))
                p.start()
                old_port = port
            time.sleep(300)


if __name__ == "__main__":
    main(sys.argv)
