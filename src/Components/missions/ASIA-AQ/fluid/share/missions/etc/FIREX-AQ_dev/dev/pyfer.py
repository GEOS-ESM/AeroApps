#! /usr/bin/env python

import sys
import yaml
import filehandler

if len(sys.argv) < 4:
    print 'Usage:', sys.argv[0], '[idate] [itime] [file(s)]'
    sys.exit(1)

cmdline = { 'idate' : int(sys.argv[1]), 'itime' : int(sys.argv[2]) }

for file in sys.argv[3:]:

    with open(file, 'r') as ymlfile:
        config = yaml.load(ymlfile)

    for protocol in config:

        fh  = filehandler.get_handler(protocol)
        cfg = config[protocol]

        for key in cfg:
    
            if not isinstance(cfg[key], dict): continue

            request = dict(cfg)
            request.update(cfg[key])
            request.update(cmdline)

            fh.get(request)
        
sys.exit(0)
