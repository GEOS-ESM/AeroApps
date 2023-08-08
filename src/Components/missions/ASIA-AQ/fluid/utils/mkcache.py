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
        if field[0] in field: return True

    return False

request = Request(interface.parse_args(sys.argv[1:]))

wx = wxservice.WXService(request)

playlist = wx.playlist()
fields   = wx.request['fields'].split(',')

for play in playlist:

    for request in play:

        if not is_requested(request['field'], fields): continue

        for r in request:

#           $field/$level/$region/$uuid.png
            rootdir = os.path.dirname(r['oname'])
            rootdir = rootdir.split('/')[0:-3] + ['.cache',r['field'],r['level'],r['region']]
            print '/'.join(rootdir)
#           lname = os.path.join(rootdir,[r['field'],r['level'],r['region'],r['uuid']])
#           print r['oname'], '--->', lname
