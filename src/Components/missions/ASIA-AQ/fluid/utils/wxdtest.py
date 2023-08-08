#! /usr/bin/env python

import os
import sys
import copy

import wxservice
import interface

from request import *

def is_requested(field, fields):

    if not fields: return True

    if len(fields) > 1:
        if field in fields: return True
    else:
        if fields[0] in field: return True

    return False

def is_excluded(field, fields):

    if not fields: return False

    for f in fields:
        if f in field: return True

    return False

request = Request(interface.parse_args(sys.argv[1:]))

wx = wxservice.WXService(request)

playlist = wx.playlist()
fields   = wx.request['fields'].split(',')
exclude  = wx.request['exclude'].split(',')

path = ['wxmapscustom','interface','field']
surface = wx.config(path+['surface','items'], [])
upper   = wx.config(path+['upper','items'], [])
fields = surface + upper

for play in playlist:

    for request in play:

        if not is_requested(request['field'], fields): continue
#       if is_excluded(request['field'], exclude): continue

        print request['region']
