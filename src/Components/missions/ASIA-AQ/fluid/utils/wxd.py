#! /usr/bin/env python

import os
import sys
import copy

import wxservice
import interface

from request import *
from taskmanager import *

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

task    = TaskManager()
request = Request(interface.parse_args(sys.argv[1:]))

wx = wxservice.WXService(request)

playlist = wx.playlist()
fields   = wx.request['fields'].split(',')
fields   = [f for f in fields if f != '']
exclude  = wx.request['exclude'].split(',')
exclude  = [f for f in exclude if f != '']

if not fields:
    path = ['wxmapscustom','interface','field']
    surface = wx.config(path+['surface','items'], [])
    upper   = wx.config(path+['upper','items'], [])
    fields = surface + upper

for play in playlist:

    for request in play:

        if not is_requested(request['field'], fields): continue
        if is_excluded(request['field'], exclude): continue

        for r in request:

#           if not os.path.isfile(r['oname']):

            cmd          = copy.copy(sys.argv[1:])
            index        = cmd.index('--start_dt')               
            start_dt     = r['time_dt']
            cmd[index+1] = start_dt.strftime('%Y%m%dT%H%M%S')
    
            cmd.append('--field '  + r['field'])
            cmd.append('--region ' + r['region'])
            cmd.append('--level '  + r['level'])
            task.spawn('wxmap.py ' + ' '.join(cmd))
            break

task.wait()
