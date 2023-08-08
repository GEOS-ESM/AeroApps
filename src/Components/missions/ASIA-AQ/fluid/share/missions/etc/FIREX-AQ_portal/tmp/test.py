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


playlist = wx.playlist()
for play in playlist:
    for req in play:
        r = next(iter(req))
        print r
