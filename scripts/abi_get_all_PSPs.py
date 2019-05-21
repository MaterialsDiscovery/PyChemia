#!/usr/bin/env python3

# This script download several sets of pseudopotentials
# form the Website of ABINIT
# The pseudopotentials will be located inside the directory
# '.abinit' of the user home directory, just aside of the
# tarballs.
# AE are the Hirshfeld All Electron (AE) densities
# This script is for Python 3 only as it depends on urllib.request
# Not backwards support for Python 2.x

import os
import ftplib
from concurrent.futures import ThreadPoolExecutor
from pychemia.code.abinit import psp_name
from urllib.request import urlopen
import tarfile
import time


def worker(filename, directory, filepath):
    ftp = ftplib.FTP('ftp.abinit.org')  # connect to host, default port
    ftp.login()  # user anonymous, passwd anonymous@
    print('Getting ftp://ftp.abinit.org/%s/%s' % (filepath, filename))
    nofile = True
    while nofile:
        try:
            ftp.retrbinary('RETR ' + filepath + filename, open(directory + '/' + filename, 'wb').write)
            nofile = False
            if os.path.getsize(directory + '/' + filename) == 0:
                os.remove(directory + '/' + filename)
        except ftplib.error_perm:
            print('Retrying download of %s' % filepath)
            ftp.close()
            time.sleep(1)
            if os.path.isfile(directory + '/' + filename):
                os.remove(directory + '/' + filename)
            ftp = ftplib.FTP('ftp.abinit.org')  # connect to host, default port
            ftp.login()  # user anonymous, passwd anonymous@
            nofile = False
        except EOFError:
            print('Retrying download of %s' % filepath)
            ftp.close()
            time.sleep(1)
            if os.path.isfile(directory + '/' + filename):
                os.remove(directory + '/' + filename)
            ftp = ftplib.FTP('ftp.abinit.org')  # connect to host, default port
            ftp.login()  # user anonymous, passwd anonymous@
            nofile = False
    ftp.close()


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
        rpath = ' https://www.abinit.org/ATOMICDATA/JTH-LDA-atomicdata.tar.gz'
    elif kind == 'PAW' and exchange == 'PBE':
        rpath = ' https://www.abinit.org/ATOMICDATA/JTH-PBE-atomicdata.tar.gz'
    elif kind == 'ONC' and exchange == 'PBE':
        rpath = ' http://www.pseudo-dojo.org/pseudos/nc-sr_pbe_standard_psp8.tgz'
    else:
        print('Not know kind of PSP')
    return rpath


def get_all_psps(basedir, exchange, kind):
    directory = basedir + os.sep + exchange + '_' + kind
    if not os.path.isdir(directory):
        os.mkdir(directory)
    if kind in ['PAW', 'ONC']:
        rpath = get_rpath_psp(kind, exchange)
        filename = rpath.split('/')[-1]
        if not os.path.isfile(directory + '/' + filename):
            u = urlopen(rpath)
            f = open(directory + os.sep + filename, 'wb')
            meta = u.info()
            file_size = int(meta.get("Content-Length"))
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
        else:
            print('All files are downloaded')

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
                print('All files are downloaded')
                break
    else:
        files = []
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
            filepath = get_rpath_psp(kind, exchange, i)
            if not os.path.isfile(directory + '/' + filename) or os.path.getsize(directory + '/' + filename) == 0:
                files.append((filename, directory, filepath))

        if len(files) == 0:
            print("All files are downloaded")
        else:
            print("Downloading %d PSP files" % len(files))
        nth = 12
        pool = ThreadPoolExecutor(nth)
        index = 0
        p = nth * [None]
        while index < len(files):
            for i in range(nth):
                if index < len(files):
                    p[i] = pool.submit(worker, *(files[index]))
                    index += 1
            for i in range(nth):
                try:
                    p[i].result()
                except AttributeError:
                    print("Complete")

        if len(files) > 0:
            print("kind == '%s' and exchange == '%s'" % (kind, exchange))
            for i in files:
                print(" %s" % str(i))


if __name__ == '__main__':

    print("Script to download PSP files from ftp.abinit.org")
    home = os.environ['HOME']
    basedir = home + "/.abinit"

    print("Files will be downloaded at: %s/.abinit" % home)

    if not os.path.isdir(basedir):
        os.mkdir(basedir)

    for exchange in ['LDA', 'GGA', 'DEN', 'PBE']:
        print('\n=> ' + exchange, end='\n')
        lista = []
        if exchange == 'LDA':
            lista = ['FHI', 'TM', 'GTH', 'PAW', 'CORE', 'HGH']
        elif exchange == 'GGA':
            lista = ['FHI', 'HGH']
        elif exchange == 'DEN':
            lista = ['AE', 'FC']
        elif exchange == 'PBE':
            lista = ['ONC', 'PAW']
        for kind in lista:
            print('--> %5s ' % kind, end=' ')
            get_all_psps(basedir, exchange, kind)
