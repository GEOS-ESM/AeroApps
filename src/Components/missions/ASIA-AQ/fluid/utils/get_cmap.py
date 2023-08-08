#! /usr/bin/env python

import os
import sys

import numpy as np

import wxservice
import interface

from request import *

def write_cmap(name, cmap):

    print '   ',name + ':'

    for key in ['cmap', 'red', 'green', 'blue', 'alpha', 'reverse', 'scale']:

        map = cmap.get(key, None)
        if not map: continue

        if isinstance(map, list):
            print ' '*6 + key + ':'

            for segment in map:
                print ' '*7, '-', segment
        else:
            print ' '*6 + key + ':', map

    print ''

request = Request(interface.parse_args(sys.argv[1:]))
wx      = wxservice.WXService(request)

colorbar = wx.config(['attribute', 'colorbar'])

for name, clist in colorbar.iteritems():

    if isinstance(clist, dict):
        write_cmap(name, clist)
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

    write_cmap(name, cmap)

sys.exit(0)
