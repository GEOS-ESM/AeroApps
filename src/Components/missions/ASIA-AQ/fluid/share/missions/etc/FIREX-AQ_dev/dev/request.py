import uuid
import datetime as dt

from string import *

novalue = object()

class Request(dict):

    def __init__(self, request, encoder=None):

        super(Request,self).__init__(request)

        self.encoder = encoder
        self.hotkeys  = ['fcst_dt', 'time_dt', 'field', 'level', 'region',
                         'geometry', 'stream', 'collection', 'basemap']

    def get_name(self, request=None, hotkeys=None, encoder=novalue):

        if request is None:
            request = self

        if hotkeys is None:
            hotkeys = self.hotkeys

        if encoder is novalue:
            encoder = self.encoder

        keystr = []
        for key in hotkeys:
            keystr.append(str(key) + ':' + repr(request.get(key,None)))

        keystr = ','.join(keystr)

        if encoder:
            keystr += ',' + encoder.encode(request)

        keystr = keystr.replace("u'", "'")
        return str(uuid.uuid3(uuid.NAMESPACE_DNS,keystr))

    def get_rname(self, request=None):

        hotkeys  = ['lat', 'lon', 'mproj', 'mpvals', 'geometry', 'basemap']

        return self.get_name(request, hotkeys)

    def get_key(self, request=None):

        if request is None:
            request = self

        out_key = ''
        for key in ['theme', 'motif', 'config']:
            if request.get(key, None):
                out_key += '+' + key + ':' + '+'.join(request[key])

        if not out_key: return 'unknown'

        return out_key

    def __iter__(self):

        t      = self['start_dt']
        end_dt = self['end_dt']
        tinc   = self['t_deltat']

        request = dict(self)

        while t <= end_dt:

            request['time_dt'] = t
            request['gtime']   = t.strftime("%Hz%d%b%Y")
            request['uuid']    = self.get_name(request)

            ftime = request['fcst_dt']
            if not ftime: ftime = t
            tau   = t - ftime
            request['tau'] = "%03d"%(int(tau.total_seconds() / 3600),)

            oname = self.get('oname','$uuid.png')
            oname = Template(oname).safe_substitute(request)
            oname = t.strftime(oname)
            oname = ftime.strftime(oname)

            request['oname'] = oname

            yield Request(request, self.encoder)

            t += tinc

    def __call__(self, key, default=None):

        value = self.get(key, default)

        if value and isinstance(value, basestring):
            return value.split('/')[-1]

        return value
