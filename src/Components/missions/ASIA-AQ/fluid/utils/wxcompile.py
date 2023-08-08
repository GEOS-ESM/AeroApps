#! /usr/bin/env python

import os
import sys
import copy
import json

import config
import wxservice
import interface

from request import *

request  = Request(interface.parse_args(sys.argv[1:]))
wx       = wxservice.WXService(request)

bin_name = request.get_key()
bin_path = wx.config('bin_path', './')
bin_path = os.path.join(bin_path, bin_name)

print bin_path

try:
    os.makedirs(bin_path, 0755)
except:
    pass

playlist = wx.playlist()

i = 0
for play in playlist:

    file = 'wx%03d.json'%(i) 
    file = os.path.join(bin_path, file)
    print file

#   wx.config['plotservice_'] = wx.ps
#   wx.config['theme_'] = wx.ps.name
    wx.config.serialize(wx.config)
    with open(file, 'w') as outfile:
        json.dump(wx.config, outfile)

    i += 1

sys.exit(0)
