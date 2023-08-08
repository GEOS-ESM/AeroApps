#! /usr/bin/env python

import os
import sys
import yaml

pathname     = os.path.dirname(sys.argv[0])
install_path = os.path.abspath(pathname)

flag = 1
if len(sys.argv[1:]) > 1: flag = int(sys.argv[2])

etc_path  = os.path.dirname(install_path) + os.sep + 'etc'
root_path = os.path.dirname(install_path)+os.sep+'share'+os.sep+'missions'
stream    = sys.argv[1]
instance  = 'FX_' + sys.argv[1]

file = root_path + os.sep + 'FIREX-AQ' + os.sep + 'aeronet.csv'
with open(file,'r') as f: lines = f.readlines()
stations = [x.strip() for x in lines if 'name,lon,lat' not in x]

file = root_path + os.sep + 'FIREX-AQ' + os.sep + 'cities.csv'
with open(file,'r') as f: lines = f.readlines()
stations += [x.strip() for x in lines if 'name,lon,lat' not in x]

file = root_path + os.sep + instance + os.sep + 'interface.yml'
with open(file,'r') as f: interface = yaml.load(f)
vars = interface['wxmapscustom']['interface']['field']['upper']['items']

file = root_path + os.sep + 'FIREX-AQ' + os.sep + 'custom.yml'
with open(file,'r') as f: config = yaml.load(f)
plots = config['wxmapscustom']['plot']

tab2 = ' '*2
tab4 = ' '*4
tab6 = ' '*6
tab8 = ' '*8

print 'wxmapscustom:'
print ''
print '  plot:'
print ''

for name in vars:

    plot = plots[name]
    NAME = name.upper()

    for i,loc in enumerate(stations):

        loc   = loc.split(',')
        lon   = float(loc[1])
        lat   = float(loc[2])
        place = loc[0].replace('.',' ')

        if flag ==  1 and lon < 0.0:   lon += 360.0
        if flag == -1 and lon > 180.0: lon -= 360.0

        if lon > 180.0:
            clon = ('%6.2fW'%(abs(lon-360.0),)).strip()
        elif lon < 0.0:
            clon = ('%6.2fW'%(abs(lon),)).strip()
        elif lon >   0.0:
            clon = ('%6.2fE'%(lon,)).strip()
        else:
            clon = ('%6.2f'%(lon,)).strip()

        if lat > 0.0:
            clat = ('%5.2fN'%(lat,)).strip()
        elif lat < 0.0:
            clat = ('%5.2fN'%(abs(lat),)).strip()
        else:
            clat = ('%5.2f'%(lat,)).strip()

        shortname = name + '_gram%03d'%(i,)

        layer = plot['layers']
        alt_layer = ['ground',name,'puffy_cloud','barbs','cloud_water','pbltop']
        if stream == 'GEOS': layer = alt_layer

        if i == 0:

            print  tab4 + shortname + ': &' + NAME
#           print  tab6 + 'long_name:', repr(plot['long_name'])
            print  tab6 + 'long_name:', repr(name)
            print  tab6 + 'header:', place + ' (' + clon + ',' + clat + ')'
            print  tab6 + 'tag:', place
            print  tab6 + "xlabel: ''"
            print  tab6 + 'levels: [0]'
            print  tab6 + 'layers: [' + ','.join(layer) + ']'
            title = plot['title']
            i = title.find('hPa')
            if i >= 0: title = title[i+4:]
            print  tab6 + "title: '" + title + "'"
            print  tab6 + 'parea:  1 10 1 7.5'
            print  tab6 + 'lev:  1000 100'
            print  tab6 + 'lat: ', lat
            print  tab6 + 'lon: ', lon
            print  tab6 + "time: '$tbeg $t5day'"
            print  tab6 + "ylab: 'Pressure (mb)'"
            print  tab6 + "shape: 'off'"

            if layer[0] in plot:
                print tab6 + layer[0] + ':'
                for k,v in plot[layer[0]].iteritems():
                    print tab8 + k + ': ' + repr(v)
            print ''

        else:

            print  tab4 + shortname + ':'
            print  tab6 + '<<: *' + NAME
#           print  tab6 + 'long_name:', plot['long_name'] + ' ' + place
            print  tab6 + 'header:', place + ' (' + clon + ',' + clat + ')'
            print  tab6 + 'tag:', place
            print  tab6 + 'lat: ', lat
            print  tab6 + 'lon: ', lon
            print ''
