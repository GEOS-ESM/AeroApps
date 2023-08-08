#! /usr/bin/env python

import sys
from taskmanager import *

if len(sys.argv) < 5:
    print 'Usage:', sys.argv[0], '[idate] [itime] [resource] [file(s)]'
    sys.exit(1)

task     = TaskManager()
template = ' '.join(sys.argv[1:4])

for file in sys.argv[4:]:
    args = template + ' ' + file
    print 'firex_ops.py ' + args
    task.spawn('firex_ops.py ' + args)

task.wait()
