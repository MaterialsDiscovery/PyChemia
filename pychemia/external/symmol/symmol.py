import os
import re
import subprocess
import numpy as np

__author__ = 'Guillermo Avendano Franco'


def get_point_group(st, executable='symmol', dcm=0.2, dcme=0.2):

    if not cmd_exists(executable):
        return 'no symmol'

    wf = open('symmol.in', 'w')
    wf.write(' %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f\n' % (1, 1, 1, 90, 90, 90))
    wf.write(' %d %d %9.6f %9.6f\n' % (1, 1, dcm, dcme))
    for i in st:
        name = i.symbols[0]
        wf.write('%6s%2d%9.5f%9.5f%9.5f\n' % tuple([name.ljust(6), 1] + list(i.position)))
    wf.close()

    rf = open('symmol.in')

    if np.max(np.abs(st.positions.flatten())) >= 100:
        return 'ERR'

    data = subprocess.check_output([executable], stdin=rf)
    data = data.decode('utf-8')
    ans = re.findall('Schoenflies symbol =([\d\w\s]+)CSM', data)
    rf.close()
    if len(ans) != 1:
        pg = 'no symmetry'
    else:
        pg = ans[0].strip()
        os.remove('symmol.in')
    return pg


def cmd_exists(cmd):
    return subprocess.call("type " + cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE) == 0
