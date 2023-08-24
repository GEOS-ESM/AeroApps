#! /usr/bin/env python

import os
import sys
import yaml
import filehandler

if len(sys.argv) < 5:
    print 'Usage:', sys.argv[0], '[idate] [itime] [resource] [file(s)]'
    sys.exit(1)

cmdline = { 'idate' : int(sys.argv[1]), 'itime' : int(sys.argv[2]) }
with open(sys.argv[3], 'r') as ymlfile: resource = yaml.load(ymlfile)

for file in sys.argv[4:]:

    with open(file, 'r') as ymlfile:
        config = yaml.load(ymlfile)

    model = os.path.basename(file)
    if 'model' not in config: config['model'] = os.path.splitext(model)[0]

    if 'shift_dt' not in cmdline:
        cmdline['shift_dt'] = int(config.get('shift_dt', 0))

    fh  = filehandler.PlotHandler(resource)

    request = dict(resource)
    request.update(config)
    request.update(cmdline)

    iret = fh.handler(request)

sys.exit(iret)
