#! /usr/bin/env python

import os
import sys
import yaml

with open('datagram.yml','r') as f: config = yaml.load(f)

stations = config['cities']
plots    = config['plot']

tab2 = ' '*2
tab4 = ' '*4
tab6 = ' '*6

for name,plot in plots.iteritems():

    NAME = name.upper()
    print '  plot:'
    print ''

    for i,loc in enumerate(stations):

        loc  = loc.split()
        lon  = float(loc[0])
        lat  = float(loc[1])
        place = (' '.join(loc[2:-1])).replace('.',' ')

        shortname = name + '%03d'%(i,)

        if i == 0:

            print  tab4 + shortname + ': &' + NAME
            print  tab6 + 'long_name:', plot['long_name']
            print  tab6 + 'header:', place
            print  tab6 + 'tag:', place
            print  tab6 + "xlabel: ''"
            levels = [str(l) for l in plot['levels']]
            print  tab6 + 'levels: [' + ','.join(levels) + ']'
            print  tab6 + 'layers: [' + ','.join(plot['layers']) + ']'
            print  tab6 + "title: '" + plot['title'] + "'"
            print  tab6 + 'parea: ', plot['parea']
            print  tab6 + 'lev: ', plot['lev']
            print  tab6 + 'lat: ', lat
            print  tab6 + 'lon: ', lon
            print  tab6 + "time: '" + plot['time'] + "'"
            print  tab6 + "ylab: '" + plot['ylab'] + "'"
            print  tab6 + "shape: '" + plot['shape'] + "'"
            print ''

        else:

            print  tab4 + shortname + ':'
            print  tab6 + '<<: *' + NAME
#           print  tab6 + 'long_name:', plot['long_name'] + ' ' + place
            print  tab6 + 'header:', place
            print  tab6 + 'tag:', place
            print  tab6 + 'lat: ', lat
            print  tab6 + 'lon: ', lon
            print ''
