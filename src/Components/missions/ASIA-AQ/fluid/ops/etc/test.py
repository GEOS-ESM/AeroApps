#! /usr/bin/env python

import sys
import gradstime as gt

if len(sys.argv)-1 != 3:
    print 'Usage:', sys.argv[0], '[idate] [itime] [string]'
    sys.exit(1)

idate = int(sys.argv[1])
itime = int(sys.argv[2])
s     = sys.argv[3]

g = gt.GradsTime(idate, itime)

print g.stritime(s, 0)
print g.strvtime(s, 0)
print g.strftime(s, 24)

sys.exit(0)
