import os
import sys
import copy
import glob
import collections
import datetime as dt

import config

from request import *
from encoder import *

class WXService(object):

    def __init__(self, request=None,
                       plotservice=None,
                       dataservice=None,
                       mapservice=None):

        self.ps      = plotservice
        self.ds      = dataservice
        self.ms      = mapservice
        self.config  = config.Config()
     
        self.counter = 0
        self.recursion_depth = 0
        self.MAX_RECURSION_DEPTH = 20

        self.register()

        if request is None:
            return

        self.request = copy.deepcopy(request)

        install_path = os.path.dirname(sys.argv[0])
        file = request.get('rc', None)
        if not file: file = os.path.join(install_path, 'wxmap.rc')

        resource = self.config.read(file)
        self.config.mount(resource, '/')

        self.configure()
        self.register(config=self.config)

        self.reset(request.get('reset',[]))
        self.configure(request.get('theme',[]), ext='.init')
        self.configure(request.get('config',[]))

        self.base_config = copy.deepcopy(self.config)

    def register(self, config=None,
                       dataservice=None,
                       plotservice=None,
                       mapservice=None):

        if config is not None:
            self.config = config
            self.base_config = copy.deepcopy(self.config)

        if dataservice is not None:
            self.ds = dataservice

        if plotservice is not None:
            self.ps = plotservice

        if mapservice is not None:
            self.ms = mapservice

        if self.config is not None and self.ds is not None:
            self.ds.register(config=self.config)

        if self.ps is not None:
            self.ps.register(config=self.config, dataservice=self.ds)

        if self.ms is not None:
            self.ms.register(config=self.config, dataservice=self.ds)

    def get_plot(self, request):

        if self.ps is None:
            return None

        return self.ps.get_plot(request)

    def get_map(self, request):

        if self.ps is None or self.ms is None:
            return None

        plot = self.ps.get_plot(request)

        return self.ms.get_maps(plot)

    def update_request(self, request, play):

        user = request.get('user',{})

        for key in play:

            is_specified = user.get(key,0)
            if is_specified: continue

            value = play[key]
            if not isinstance(value, basestring): continue

            try:
                request[key] = dt.datetime.strptime(value,'%Y%m%dT%H%M%S')
            except:
                request[key] = value

            if key == 'time_dt':
                request['start_dt'] = request[key]
                request['end_dt']   = request[key]

    def configure(self, cfg=None, ext='.yml'):

        self.recursion_depth += 1

        assert self.recursion_depth <= self.MAX_RECURSION_DEPTH, \
               'Too many nested configuration references or circular reference'

        if cfg is None:
            names = self.config.get('config',[])
        elif isinstance(cfg, list):
            names = cfg
        elif isinstance(cfg, dict):
            return
        else:
            names = [cfg]

        paths  = self.config.get('config_path',None)

        if not paths:
            paths = [os.getcwd()]
        else:
            paths = paths.split(':') + [os.getcwd()]

        for name in names:

            name, srch_path = self.provenance(name, ext, paths)
                
            for path in srch_path:
                for file in self.list(os.path.join(path,name),ext):
                    resource = self.config.read(file)

                    self.reset(resource.get('reset',[]))
                    self.configure(resource.get('theme',[]),ext='.init')
                    self.configure(resource.get('config',[]))

                    self.config.mount(resource)
                    ps = self.get_plotservice(resource)
                    if ps: self.register(plotservice=ps)

        self.recursion_depth -= 1

    def provenance(self, name, ext, paths):

        file = name
        if ext not in name: file += ext

        path = os.path.dirname(name)
        if not path: path = os.getcwd()

        if os.path.isdir(name) or os.path.isfile(file):
            return (os.path.basename(name), [path])
        else:
            return (name, paths)

    def get_plotservice(self, resource):

        for key in resource.keys():
            if isinstance(resource[key],dict):
                service = resource[key].get('service',None)
                if getattr(service, 'get_plot', None):
                    service.set_name(key)
                    return service

        return None

    def list(self, node, ext='.yml'):

        files = []

        f = node
        if ext not in node: f += ext

        if os.path.isfile(f):
            files.append(f)

        if os.path.isdir(node):

            listing = glob.glob(os.path.join(node,'*'+ext))
            files  += [f for f in listing if os.path.isfile(f)]

        return files

    def reset(self, reset):

        for name in reset:
            self.config[name] = {}

    def playlist(self):

        playlist = self.base_config.get('playlist',None)
        if not playlist: playlist = {'default':{}}

        for name in playlist:

            play   = playlist[name]
            config = copy.deepcopy(self.base_config)

            self.register(config=config)
            self.configure(play.get('theme',[]),ext='.init')
            self.configure(play.get('config',[]))
            self.configure(self.request.get('config',[]))

            self.config['plotservice_'] = self.ps
            self.config['theme_'] = self.ps.name
            self.config['play_']  = play

            yield self.play(play)

    def play(self, play):

        response = copy.deepcopy(self.request)
        self.update_request(response, play)

        r_field  = response.get('field', None)
        r_region = response.get('region', None)
        r_level  = response.get('level', None)

        theme   = self.ps.name
        fields  = self.config([theme, 'plot'], {}).keys()
        regions = self.config(['region'], {}).keys()

        rmask  = self.config('regions', regions)
        rmask  = play.get('regions', rmask)

        fmask  = self.config('fields', fields)
        fmask  = play.get('fields', fmask)

        for field in fields:

            if r_field and r_field != field: continue
            if field not in fmask: continue

            response['field'] = field
            long_name  = self.config([theme, 'plot', field, 'long_name'], '')
            tag_name   = self.config([theme, 'plot', field, 'tag'], '')
            levels     = self.config([theme, 'plot', field, 'levels'],[0])
            levels     = self.config('levels', levels)

            lmask  = self.config('levels', levels)
            lmask  = play.get('levels', lmask)

            response['field_name'] = '_'.join(long_name.split())
            response['tag_name'] = '_'.join(tag_name.split())

            for region in regions:

                if r_region and r_region != region: continue
                if region not in rmask: continue
                response['region'] = region

                long_name = self.config(['region', region, 'long_name'], '')
                response['region_name'] = '_'.join(long_name.split())

                for level in levels:

                    if r_level and r_level != str(level): continue
                    if level not in lmask: continue
                    response['level']  = str(level)

                    yield Request(response, Encoder(self.ps))

    def get_capabilities(self):

        streams = {}
        fields  = {}
        regions = {}

        for p in self.playlist():

            play = self.config['play_']

            root    = [self.ps.name]

            path    = ['stream']
            result  = self.config.find(path, name='long_name',depth=1)
            result  = self.config.get_items(result)
            default = result.keys()
            default = self.config('streams', default)
            mask    = play.get('streams', default)
            streams.update({ k:v for k,v in result.iteritems() if k in mask })

            path      = root + ['plot']
            result    = self.config.find(path, name='long_name',depth=1)
            result    = self.config.get_items(result)
            f_default = result.keys()
            f_default = self.config('fields', f_default)
            mask      = play.get('fields', f_default)
            fields.update({ k:v for k,v in result.iteritems() if k in mask })
            
            path    = ['region']
            result  = self.config.find(path, name='long_name',depth=1)
            result  = self.config.get_items(result)
            default = result.keys()
            default = self.config('regions', default)
            mask    = play.get('regions', default)
            regions.update({ k:v for k,v in result.iteritems() if k in mask })
            
#           default = fields.keys()[0]
            default = f_default[0]
            field   = self.request.get('field',default)
            if not field: field = default

            path    = root + ['plot',field]
            result  = self.config.find(path, name='levels',depth=1)
            result  = self.config.get_values(result)[0]
            default = result
            default = self.config('levels', default)
            mask    = play.get('levels', default)
            levels  = [lev for lev in result if lev in mask]
            
        return dict(zip(['stream', 'field', 'region', 'level'],
                   [streams,  fields,  regions,  levels]))

    def get_user_interface(self):

        response  = self.get_capabilities()
        interface = self.config([self.ps.name, 'interface'], {})

        ui = collections.OrderedDict()

        for section_name in ['stream', 'field', 'region', 'level']:

            section    = interface.get(section_name, {})
            od_section = collections.OrderedDict()
            default    = response.get(section_name, {})

            ui[section_name] = od_section

            for group_name in section.get('groups', []):

                group = section.get(group_name, {})
                od_group = collections.OrderedDict()
                od_item  = collections.OrderedDict()
                od_group.update(group)

                for item in group.get('items', []):

                    if isinstance(item, dict):
                        od_item.update(item)
                    else:
                        od_item[item] = default.get(item, 'Unknown')

                od_group['items'] = od_item
                od_section[group_name] = od_group

            if not od_section:
                od_group = collections.OrderedDict()
                od_group['items'] = default
                od_section['default'] = od_group

        return ui

    def renew(self, freq=None):

        if freq is None or freq <= 0:
            return self.ds

        if self.ds is None:
            return self.ds

        self.counter += 1
        if (self.counter % freq) == 0:
            ds = type(self.ds)()
            self.register(dataservice = ds)
            self.counter = 0

        return self.ds

    def clone(self):

        wx = WXService()

        if self.config: wx.config = copy.deepcopy(self.config)
        if self.base_config: wx.base_config = copy.deepcopy(self.base_config)
        if self.request: wx.request = copy.deepcopy(self.request)

        if self.ps: wx.ps = type(self.ps)(name=self.ps.name)
        if self.ds: wx.ds = type(self.ds)()
        if self.ms: wx.ms = type(self.ms)()

        wx.register()

        return wx

class WXServiceLite(WXService):

    def __init__(self, request=None,
                       plotservice=None,
                       dataservice=None,
                       mapservice=None):

        self.ps      = plotservice
        self.ds      = dataservice
        self.ms      = mapservice
        self.config  = config.Config()
     
        self.counter  = 0
        self.num_play = 0

        self.register()

        if request is None:
            return

        r = copy.deepcopy(request)
        self.request = Request(r)

        install_path = os.path.dirname(sys.argv[0])
        file = self.request.get('rc', None)
        if not file: file = os.path.join(install_path, 'wxmap.rc')
        resource = self.config.read(file)

        self.bin_path  = resource.get('bin_path', './')
        self.bin_path  = os.path.join(self.bin_path, self.request.get_key())

        listing = glob.glob(os.path.join(self.bin_path,'*.json'))
        self.num_play = len(listing)

        base_file = os.path.join(self.bin_path, 'wx000.json')

        cfg = config.Config(self.config.readJSON(base_file))
        plotservice = cfg['plotservice_']
        plotservice.set_name(cfg['theme_'])
        self.register(config=cfg, plotservice=plotservice)

    def playlist(self):

        for i in range(0,self.num_play):

            file   = 'wx%03d.json'%(i)
            file   = os.path.join(self.bin_path, file)
            cfg    = config.Config(self.config.readJSON(file))

            plotservice = cfg['plotservice_']
            plotservice.set_name(cfg['theme_'])
            play = cfg['play_']
            self.register(config=cfg, plotservice=plotservice)

            yield self.play(play)
