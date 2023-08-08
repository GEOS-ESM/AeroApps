#! /usr/bin/env python

import os
import sys
import yaml

model = sys.argv[1]

pathname     = os.path.dirname(sys.argv[0])
install_path = os.path.abspath(pathname)

etc_path  = os.path.dirname(install_path) + os.sep + 'etc'
root_path = os.path.dirname(install_path)+os.sep+'share'+os.sep+'missions'
instance  = sys.argv[1]

# Read in fields for model

file = root_path + os.sep + 'FIREX-AQ' + os.sep + 'field.yml'
with open(file,'r') as f: config = yaml.load(f)
fields = config['field'][model]
fields = [k[1:] for k in fields.keys()]

# Read in plots

file = root_path + os.sep + 'FIREX-AQ' + os.sep + 'custom.yml'
with open(file,'r') as f: config = yaml.load(f)
plot = config['wxmapscustom']['plot']
plot_names  = plot.keys()

levels = {}
for name in plot_names:
    p = plot[name]
    levels[name] = p['levels']

surface = []
upper = []

for name in fields:

    if name not in plot_names: continue

    if len(levels[name]) == 1:
        surface.append(name)
    else:
        upper.append(name)

tab2  = ' ' * 2
tab4  = ' ' * 4
tab6  = ' ' * 6
tab8  = ' ' * 8
tab10 = ' ' * 10

print 'wxmapscustom:'
print ''
print tab2 + 'interface:'
print ''
print tab4 + 'field:'
print ''
print tab6 + 'long_name: VARIABLES'

groups = []
if surface: groups += ['surface']
if upper: groups += ['upper']

print tab6 + 'groups: ' + repr(groups)
print ''

if surface:

    print tab6 + 'surface:'
    print tab8 + 'long_name: Surface'
    print tab8 + 'type: button'
    print tab8 + 'items:'
    for item in sorted(surface): print tab10 + '- ' + item
    print ''

if upper:

    print tab6 + 'upper:'
    print tab8 + 'long_name: Pressure Level'
    print tab8 + 'type: button'
    print tab8 + 'items:'
    for item in sorted(upper): print tab10 + '- ' + item
