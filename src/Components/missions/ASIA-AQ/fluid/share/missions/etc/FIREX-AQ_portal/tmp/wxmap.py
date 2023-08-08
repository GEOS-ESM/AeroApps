#! /usr/bin/env python

import os
import sys
import copy
import wxservice
import interface
import gradsdataservice as dataservice
import gradsmapservice as mapservice

from request import *

request = Request(interface.parse_args(sys.argv[1:]))

#wx = wxservice.WXService(request)
wx = wxservice.WXServiceLite(request)

ds = dataservice.Service()
ms = mapservice.Service()

wx.register(dataservice = ds)
wx.register(mapservice  = ms)

playlist = wx.playlist()

for play in playlist:

    for request in play:

        for r in request:
            opath = os.path.dirname(r['oname'])
            try:
                os.makedirs(opath, 0755)
            except:
                pass

            ds  = wx.renew(10)
            map = wx.get_map(r)
            print map
