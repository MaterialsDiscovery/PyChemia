#!/usr/bin/env python

import os
import subprocess
import pychemia

for ifile in [x for x in os.listdir('.') if x[-5:] == 'ascii']:
    st = pychemia.io.ascii.load(ifile)
    pov = pychemia.gui.StructurePovray(st)
    povfile = ifile[:-5] + 'pov'
    jpgfile = ifile[:-5] + 'jpg'
    pov.write_povray(povfile)
    subprocess.call(['povray', '+I%s' % povfile, '+FJ', '+O%s' % jpgfile, '+W1600', '+H1200', '+Q9'])
