#! /usr/bin/env python

import os
import sys
import yaml

pathname     = os.path.dirname(sys.argv[0])
install_path = os.path.abspath(pathname)

gflag = 1
uflag = 1
if len(sys.argv[1:]) > 1:

    uflag = int(sys.argv[2])
    gflag = int(sys.argv[-1])

etc_path  = os.path.dirname(install_path) + os.sep + 'etc'
root_path = os.path.dirname(install_path)+os.sep+'share'+os.sep+'missions'
stream    = sys.argv[1]
instance  = 'FX_' + stream

stations = {}

file = root_path + os.sep + 'FIREX-AQ' + os.sep + 'aeronet.csv'
with open(file,'r') as f: lines = f.readlines()
stations['aeronet'] = [x.strip() for x in lines if 'name,lon,lat' not in x]

file = root_path + os.sep + 'FIREX-AQ' + os.sep + 'cities.csv'
with open(file,'r') as f: lines = f.readlines()
stations['cities'] = [x.strip() for x in lines if 'name,lon,lat' not in x]

#file = root_path + os.sep + instance + os.sep + 'region.yml'
#with open(file,'r') as f: interface = yaml.load(f)
#vars = interface['wxmapscustom']['interface']['field']['upper']['items']

#file = root_path + os.sep + 'FIREX-AQ' + os.sep + 'region.yml'
#with open(file,'r') as f: config = yaml.load(f)
#plots = config['wxmapscustom']['plot']

tab2 = ' '*2
tab4 = ' '*4
tab6 = ' '*6

URL = 'https://portal.nccs.nasa.gov/datashare/gmao_ops/pub/fp/.internal/FIREX-AQ/$fcst/$stream/datagrams/$field/model.$stream.$fcst.000.$field.0'

url = URL.replace('$stream', stream)
if gflag == 0: stations = {}

for type, locations in stations.iteritems():

    TYPE = type.upper()
    print type + ':' + ' &' + TYPE 
    print ''

    for i,loc in enumerate(locations):

        loc   = loc.split(',')
        lon   = float(loc[1])
        lat   = float(loc[2])
        place = loc[0].replace('.',' ')
        _place = place.replace(' ','_')

        if gflag ==  1 and lon < 0.0:   lon += 360.0
        if gflag == -1 and lon > 180.0: lon -= 360.0

        if uflag == 0:
            print '  - ', lon, lat
        else:
            print '  - ', lon, lat, place, url + '.' + _place + '.png'

    print ''

file = root_path + os.sep + 'FIREX-AQ' + os.sep + 'region.yml'
with open(file,'r') as f: lines = f.readlines()
lines = [line.rstrip() for line in lines]

for line in lines:
    print line
    if 'addlayers' in line and stations:
        print ''
        print '    aeronet: *AERONET'
        print '    cities: *CITIES'
