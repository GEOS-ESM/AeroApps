#! /usr/bin/env python

import sys
import gradsdataservice as dataservice

ds = dataservice.Service()

try:
    ret = ds.open(sys.argv[1])
except:
    print 'failed'
