#! /usr/bin/env python

import sys
import yaml
import filehandler
import datetime as dt

if len(sys.argv) < 5:
    print 'Usage:', sys.argv[0], '[idate] [itime] [resource] [file(s)]'
    sys.exit(1)

cmdline = { 'idate' : int(sys.argv[1]), 'itime' : int(sys.argv[2]) }
with open(sys.argv[3], 'r') as ymlfile: resource = yaml.load(ymlfile)

idate   = int(sys.argv[1])
current = dt.datetime.utcnow()
cdate   = int(current.strftime('%Y%m%d'))

for file in sys.argv[4:]:

    status = 0

    with open(file, 'r') as ymlfile:
        config = yaml.load(ymlfile)

    if 'shift_dt' not in cmdline:
        cmdline['shift_dt'] = int(config.get('shift_dt', 0))

    handlers = config.get('handlers',None)
    if handlers == '': continue

    if handlers:
        handlers = handlers.split(':')
    else:
        handlers = [k for k,v in config.iteritems() if isinstance(v,dict)]

    for name in handlers:

        fh  = filehandler.get_handler(name, resource, status)
        if not fh: continue

        cfg    = config[name]
        groups = [v for k,v in cfg.iteritems() if isinstance(v,dict)]
        if not groups: groups = [cfg]

        for group in groups:
    
            request = dict(cfg)
            request.update(group)
            request.update(cmdline)

            iret = fh.handler(request)
            if iret != 0: status = 2

        print name, ': status = ', status

if idate < cdate: status = 0
sys.exit(status)
