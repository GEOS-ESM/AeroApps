import os
import sys
import fnmatch
from string import *
import datetime as dt

class DataService(object):

    def __init__(self, config=None):

        self.config  = config
        self.lats    = None
        self.lons    = None

    def get_capabilities(self, request):

        fh = self.open(request)

        info = {}
        info['time']  = fh.time
        info['field'] = fh.field

        return info

    def register(self, config=None):

        if config is not None:
            self.config = config

    def inquire(self, request, **kwargs):

        vtime = request['time_dt']
        ftime = request.get('fcst_dt',None)

        if ftime is None: ftime = vtime

        stream = request['stream']
        uri = self.config['stream'][stream]['uri']
        uri = uri.split('(')[0]
        uri = vtime.strftime(uri)
        uri = ftime.strftime(uri)
      
        return Template(uri).safe_substitute(request)

    def pad_region(self, pad=0):

        qf = self.query("file")
        qd = self.query("dims")

#       print 'Dimension Information'
#       print qd.dfile
#       print qd.x, qd.nx, qf.nx
#       print qd.y, qd.ny, qf.ny

        if pad > 0:

          x1 = int(qd.x[0] - pad)
          x2 = int(qd.x[1] + pad)
          y1 = int(qd.y[0] - pad)
          y2 = int(qd.y[1] + pad)

          if y1 < 0: y1 = 1
          if y2 > qf.ny: y2 = qf.ny

          if (x2-x1+1) <= qf.nx:
              self.cmd("set x %d %d"%(x1, x2))
#             print "set x %d %d"%(x1, x2)

          self.cmd("set y %d %d"%(y1, y2))
#         print "set y %d %d"%(y1, y2)

          self.lats = qd.y
          self.lons = qd.x

    def reset_region(self):

        if self.lats:
            self.cmd("set y %d %d"%self.lats)
#           print "set y %d %d"%self.lats
            self.lats = None

        if self.lons:
            self.cmd("set x %d %d"%self.lons)
#           print "set x %d %d"%self.lons
            self.lons = None

class FileHandle(object):

    def __init__(self, file, fileinfo, diminfo, timeinfo, ctlinfo, **kwargs):

        self.__dict__.update(kwargs)

        self.file     = file
        self.fileinfo = fileinfo
        self.timeinfo = timeinfo
        self.diminfo  = diminfo
        self.ctlinfo  = ctlinfo

        self.time  = {}
        self.field = {}

        self.time['nt']       = fileinfo.nt
        self.time['start_dt'] = timeinfo.tyme1
        self.time['end_dt']   = timeinfo.tyme2
        self.time['time_dt']  = timeinfo.tyme1

        for varname, nlevels, dimnames, vartitle in fileinfo.var_info:

            var     = {}
            varname = varname.lower()

            var['long_name'] = vartitle
            var['level']     = ctlinfo.zlevs
            var['units']     = ''

            self.field[varname] = var

class GEOSDDF(object):

    def __init__(self):
        pass

    def query(self, config, fmt="%Y%m%d_%H"):

        collections = {}
        times       = {}

        name = os.path.basename(config['uri'])
        path = os.path.dirname(config['uri'])
        root = os.path.dirname(path)

        for (root, dirs, files) in os.walk(root):

            for f in files:

                nodes = f.split('.')

                if len(nodes) != 2:
                    continue

                collection = nodes[0]
                date_time  = nodes[1]

                try:
                    time_dt = dt.datetime.strptime(date_time, fmt)
                except:
                    continue

                time_str = time_dt.strftime("%Y%m%dT%H%M%S")

                collections[collection] = 1
                times[time_str] = time_dt

        collections = sorted(collections.keys())
        times       = sorted(times.keys())

        return (times, collections)
