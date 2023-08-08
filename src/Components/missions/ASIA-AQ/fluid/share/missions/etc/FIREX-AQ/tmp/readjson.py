#! /usr/bin/env python

import os
import sys
import json

if len(sys.argv)-1 != 1:
    print 'Usage:', sys.argv[0], '[json file]'
    sys.exit(1)

record = {}
file = sys.argv[1]

with open(file,'r') as f:
    for line in f:
        station = json.loads(line)

        for k,v in station.iteritems():
            k = k.decode("ascii", "ignore")
            station[k] = v

        if 'coordinates' in station:
            lat = station['coordinates']['latitude']
            lon = station['coordinates']['longitude']
            if lon < 0.0: lon += 360.0
            key = '(' + str(lon) + ' ' + str(lat) + ')'

            try:
                city = station['city'].decode("ascii", "ignore")
            except:
                city = key

            key = str(lon) + ',' + str(lat)
            if key in record: continue

            record[key] = city

for k,v in record.iteritems(): print k + ',' + v
