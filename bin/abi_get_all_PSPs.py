#!/usr/bin/env python

# This script download several sets of pseudopotentials
# form the Website of ABINIT
# The pseudopotentials will be located inside the directory
# '.abinit' of the user home directory, just aside of the
# tarballs.
# AE are the Hirshfeld All Electron (AE) densities

from __future__ import print_function
import os
import ftplib

try:
    from urllib2 import urlopen
except ImportError:
    from urllib.request import urlopen
import tarfile
import time

from pychemia import HAS_SCIPY

if HAS_SCIPY:
    from pychemia.code.abinit import psp_name
else:
    raise ImportError('scipy could not be found')


def get_rpath_psp(kind, exchange, atomicnumber=None):
    """
    Get the ftp path to get the PSP

    :param kind: (str) Source of Pseudopotentials
    :param exchange: (str) 'LDA' or 'GGA'
    :param atomicnumber: (int) Atomic Number
    :return:
    """
    rpath = None
    if kind == 'FHI' and exchange == 'LDA':
        rpath = '/pub/abinitio/Psps/LDA_FHI/'
    elif kind == 'FHI' and exchange == 'GGA':
        rpath = '/pub/abinitio/Psps/GGA_FHI/'
    elif kind == 'GTH' and exchange == 'LDA':
        rpath = '/pub/abinitio/Psps/LDA_GTH/'
    elif kind == 'CORE' and exchange == 'LDA':
        rpath = '/pub/abinitio/Psps/LDA_Core/'
    elif kind == 'HGH' and exchange == 'LDA':
        rpath = '/pub/abinitio/Psps/LDA_HGH/'
    elif kind == 'HGH' and exchange == 'GGA':
        rpath = '/pub/abinitio/Psps/GGA_HGH/'
    elif kind == 'TM' and exchange == 'LDA':
        if atomicnumber is None:
            raise ValueError('Atomic number required')
        rpath = '/pub/abinitio/Psps/LDA_TM.psps/' + str(atomicnumber).zfill(2) + '/'
    elif kind == 'AE' and exchange == 'DEN':
        rpath = '/pub/abinitio/Psps/AE_DEN/'
    elif kind == 'FC' and exchange == 'DEN':
        rpath = '/pub/abinitio/Psps/FC_DEN/'
    elif kind == 'PAW' and exchange == 'LDA':
        rpath = 'http://www.abinit.org/downloads/PAW2/ATOMICDATA/JTH-LDA-atomicdata.tar.gz'
    elif kind == 'PAW' and exchange == 'GGA':
        rpath = 'http://www.abinit.org/downloads/PAW2/ATOMICDATA/JTH-PBE-atomicdata.tar.gz'
    else:
        print('Not know kind of PSP')
    return rpath


def get_all_psps(basedir, exchange, kind):
    directory = basedir + os.sep + exchange + '_' + kind
    if not os.path.isdir(directory):
        os.mkdir(directory)
    if kind == 'PAW':
        rpath = get_rpath_psp(kind, exchange)
        filename = rpath.split('/')[-1]
        if not os.path.isfile(directory + '/' + filename):
            u = urlopen(rpath)
            f = open(directory + os.sep + filename, 'wb')
            meta = u.info()
            file_size = int(meta.get("Content-Length")[0])
            print("Downloading: %s Bytes: %s" % (filename, file_size))
            file_size_dl = 0
            block_sz = 8192
            while True:
                readed_buffer = u.read(block_sz)
                if not readed_buffer:
                    break
                file_size_dl += len(readed_buffer)
                f.write(readed_buffer)
                status = r"%10d  [%3.2f%%]" % (file_size_dl, file_size_dl * 100. / file_size)
                status += chr(8) * (len(status) + 1)
                print(status, end='')
            f.close()
            print('\n')
        try:
            tar = tarfile.open(directory + '/' + filename, 'r:gz')
            for item in tar:
                if not os.path.exists(item.name):
                    tar.extract(item, path=directory)
        except tarfile.ReadError:
            name = os.path.basename(filename)
            print(name[:name.rfind('.')], '<filename>')

    elif kind == 'HGH':
        while True:
            succeed = True
            ftp = ftplib.FTP('ftp.abinit.org')  # connect to host, default port
            ftp.login()  # user anonymous, passwd anonymous@
            ftp.cwd('pub/abinitio/Psps/LDA_HGH/')
            for filename in ftp.nlst():
                if not os.path.exists(directory + '/' + filename) or os.path.getsize(directory + '/' + filename) == 0:
                    print('Getting %s' % filename)
                    try:
                        ftp.retrbinary('RETR ' + filename, open(directory + '/' + filename, 'wb').write)
                    except ftplib.error_perm:
                        print('Failed to get ', filename)
                        succeed = False
            ftp.close()
            if succeed:
                break
    else:
        ftp = ftplib.FTP('ftp.abinit.org')  # connect to host, default port
        ftp.login()  # user anonymous, passwd anonymous@
        missing_psps = []
        for i in range(1, 113):
            if kind == 'GTH' and i > 17:
                continue
            if kind == 'CORE' and i not in [6, 7]:
                continue
            if kind == 'FHI' and exchange == 'LDA' and i in [57, 59, 63, 65, 66, 67, 90, 94, 102, 110, 111, 112]:
                continue
            if kind == 'TM' and exchange == 'LDA' and i in [104, 105, 106, 107, 108, 109, 110, 111, 112]:
                continue
            if kind == 'FHI' and exchange == 'GGA' and i in [57, 59, 63, 65, 66, 67, 90, 94, 102, 110, 111, 112]:
                continue
            if kind == 'AE' and exchange == 'DEN' and i in [57, 59, 63, 65, 66, 67, 90, 94, 102, 110, 111, 112]:
                continue
            if kind == 'GTH' and exchange == 'LDA' and i in [2, 10]:
                continue
            if kind == 'FC' and exchange == 'DEN' and i in [63, 65, 66, 67, 110, 111, 112]:
                continue
            filename = psp_name(i, exchange, kind)
            if not os.path.isfile(directory + '/' + filename) or os.path.getsize(directory + '/' + filename) == 0:
                print('Getting...' + filename)
                nofile = True
                while nofile:
                    retr = 'RETR ' + get_rpath_psp(kind, exchange, i) + filename
                    try:
                        ftp.retrbinary(retr, open(directory + '/' + filename, 'wb').write)
                        nofile = False
                        if os.path.getsize(directory + '/' + filename) == 0:
                            os.remove(directory + '/' + filename)
                    except ftplib.error_perm:
                        print('Could not download ' + retr)
                        missing_psps.append(i)
                        ftp.close()
                        time.sleep(5)
                        print('Reconnecting...')
                        if os.path.isfile(directory + '/' + filename):
                            os.remove(directory + '/' + filename)
                        ftp = ftplib.FTP('ftp.abinit.org')  # connect to host, default port
                        ftp.login()  # user anonymous, passwd anonymous@
                        nofile = False
        ftp.close()
        if len(missing_psps) > 0:
            print("kind == '%s' and exchange == '%s' and i in %s" % (kind, exchange, missing_psps))


if __name__ == '__main__':

    home = os.environ['HOME']
    basedir = home + "/.abinit"

    if not os.path.isdir(basedir):
        os.mkdir(basedir)

    for exchange in ['LDA', 'GGA', 'DEN']:
        print('=> ' + exchange)
        lista = []
        if exchange == 'LDA':
            lista = ['FHI', 'TM', 'GTH', 'PAW', 'CORE', 'HGH']
        elif exchange == 'GGA':
            lista = ['FHI', 'HGH', 'PAW']
        elif exchange == 'DEN':
            lista = ['AE', 'FC']

        for kind in lista:
            print('--> ' + kind)
            get_all_psps(basedir, exchange, kind)
