import os
import re
import subprocess
import numpy as np

__author__ = 'Guillermo Avendano Franco'


def get_point_group(st):
    wf = open('symmol.in', 'w')
    wf.write(' %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f\n' % (1, 1, 1, 90, 90, 90))
    wf.write(' %d %d %9.6f %9.6f\n' % (1, 1, 0.2, 0.2))
    for i in st:
        name = i.symbols[0]
        wf.write('%6s%2d%9.5f%9.5f%9.5f\n' % tuple([name.ljust(6), 1] + list(i.position)))
    wf.close()

    rf = open('symmol.in')

    if np.max(np.abs(st.positions.flatten())) >= 100:
        return 'ERR'

    data = subprocess.check_output(['./symmol.exe'], stdin=rf)
    ans = re.findall('Schoenflies symbol =([\d\w\s]+)CSM', data)
    rf.close()
    if len(ans) != 1:
        pg = 'x'
    else:
        pg = ans[0].strip()
        os.remove('symmol.in')
    return pg
