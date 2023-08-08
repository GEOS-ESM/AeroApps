#! /usr/bin/env python

import os
import sys
import yaml
import uuid

import numpy as np

def write_cmap(file, cmap):

    cm_addr   = {}
    cm_ref    = {}
    cm_anchor = {}
    cm_keys   = ['cmap', 'red', 'green', 'blue', 'alpha', 'reverse', 'scale']

    for name in sorted(cmap):

        cm = cmap[name]
        id = str(uuid.uuid3(uuid.NAMESPACE_DNS,str(cm)))

        if id in cm_addr:
            cm_ref[name] = cm_addr[id]
            cm_anchor[cm_addr[id]] = 1
        else:
            cm_addr[id] = name

    print 'attribute:'
    print ''
    print '  colorbar:'
    print ''

    for name in sorted(cmap):

        if name in cm_ref:
            print '   ',name + ': *ref_' + cm_ref[name]
            print ''
            continue
        elif name in cm_anchor:
            print '   ',name + ': &ref_' + name
        else:
            print '   ',name + ':'

        cm = cmap[name]

        for key in cm_keys:

            map = cm.get(key, None)
            if not map: continue

            if isinstance(map, list):
                print ' '*6 + key + ':'

                for segment in map:
                    print ' '*7, '-', segment
            else:
                print ' '*6 + key + ':', map

        print ''

def write_rgba(file, cmap):

    cm_addr   = {}
    cm_ref    = {}
    cm_anchor = {}

    for name in sorted(cmap):

        cm = cmap[name]
        id = str(uuid.uuid3(uuid.NAMESPACE_DNS,str(cm)))

        if id in cm_addr:
            cm_ref[name] = cm_addr[id]
            cm_anchor[cm_addr[id]] = 1
        else:
            cm_addr[id] = name

    print 'attribute:'
    print ''
    print '  colorbar:'
    print ''

    for name in sorted(cmap):

        if name in cm_ref:
            print '   ',name + ': *ref_' + cm_ref[name]
            print ''
            continue
        elif name in cm_anchor:
            print '   ',name + ': &ref_' + name
        else:
            print '   ',name + ':'

        cm     = cmap[name]
        red    = cm['red']
        green  = cm['green']
        blue   = cm['blue']
        alpha  = cm['alpha']

        for index in range(len(red)):

            r = [float(c)*255.0 for c in red[index].split() if c != ' ']
            g = [float(c)*255.0 for c in green[index].split() if c != ' ']
            b = [float(c)*255.0 for c in blue[index].split() if c != ' ']
            a = [float(c)*255.0 for c in alpha[index].split() if c != ' ']

            print ' '*7, '-', "%5.3f %5.3f %5.3f %5.3f"%(r[1],g[1],b[1],a[1]) 

        print ''

########
# Main #
########

for file in sys.argv[1:]:

    CMAP = {}

    with open(file, 'r') as ymlfile:
        config = yaml.load(ymlfile)

    colorbar = config['attribute']['colorbar']

    for name, clist in colorbar.iteritems():

        if isinstance(clist, dict):
            CMAP[name] = clist
            continue

        cmap   = {}
        colors = []

        for color in clist:

            rgba  = [ float(c)/255.0 for c in color.split() if c != ' ' ]
            if len(rgba) < 4: rgba.append(1.0)
            colors.append(rgba)

        data = np.linspace(0.0, 1.0, len(colors))

        for i,channel in enumerate(['red', 'green', 'blue', 'alpha']):

            cmap[channel] = []

            for index, rgba in enumerate(colors):
                x  = data[index]
                y0 = rgba[i]
                y1 = y0
                values = '%5.3f'%(x) + ' ' + '%5.3f'%(y0) + ' ' + '%5.3f'%(y1)
                cmap[channel].append(values)

        CMAP[name] = cmap

    write_cmap(file, CMAP)

sys.exit(0)
