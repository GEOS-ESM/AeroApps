import re
import os
import copy
import sys
import json
import math
import collections
from string import *

import numpy as np

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mc

import gdsvml
import field
import evaluator
import toolkit

from mygrads.gacm import *

novalue = object()

class Plot(object):

    def __init__(self, request, service=None, passive=False, **kwargs):

        self.ds       = None
        self.config   = None
        self.theme    = None
        self.service  = None
        self.options  = dict(kwargs)

        if service:
            self.ds      = service.ds
            self.config  = service.config
            self.theme   = service.name
            self.service = service

        if passive: self.ds = None

        self.request = dict(request)
        self.fields  = []
        self.stack   = {}
        self.cmds    = None
        self.handle  = PlotHandle(**self.request)
        self.defs    = self.handle.__dict__
        self.lang    = gdsvml.GDSVML()
        self.color   = None
        self.eval    = evaluator.Evaluator(self.request, self.config, self.ds)
        self.special = { 'fileID': '$#', 'color': '$*'}
        self.lats    = None
        self.lons    = None
        self.passive = passive
        self.counter = 0
        self.panel   = request.get('panel', 0)

    def cmd(self, commands, layer=None, **kwargs):

        defs = dict(self.request)
        defs.update(self.defs)
        defs.update(kwargs)
        defs = { k:str(v) for k,v in defs.iteritems() }

        cmds = self.rsubstitute(commands, **defs)
        cmds = cmds.strip().split('\n')
        cmds = [re.sub("\s+"," ",cmd.strip()) for cmd in cmds if len(cmd) > 0]
        cmds = self.filter(cmds, self.defs)

        layer       = int(kwargs.get('zorder',0))
        self.cmds   = self.get_layer_stack('cmds', layer)
        layer_names = self.get_layer_stack('layer_names', layer)

        layer_name  = kwargs.get('layer_name', 'unknown')

        for cmd in cmds:

            directive = cmd.split(' ')[0]

            if self.lang.is_display(cmd):
                
                expr       = []
                expression = ''.join(cmd.split(' ')[1:])

                for component in expression.split(';'):
                    fld  = field.Field(component)
                    fld  = self.eval.evaluate(fld)
                    expr.append(fld)
                    layer_names.append(layer_name)

                self.cmds.append('display ' + ';'.join(expr))

                c = self.lang.get_color('set rgb $* 0 0 0 0')
                self.cmds.append('set rgb ' + str(c) + ' 0 0 0 0')
                self.cmds.append('set xlopts ' + str(c) + ' 1')
                self.cmds.append('set ylopts ' + str(c) + ' 1')

                if self.passive: continue
                self.substitute(fileID=self.ds.fileID)

            elif self.lang.is_define(cmd):

                assignment = ''.join(cmd.split(' ')[1:])
                index      = assignment.index('=')
                var        = assignment[0:index]
                expr       = assignment[index+1:]

                expr = self.eval.evaluate(field.Field(expr))

                self.eval.define(var)

                self.cmds.append('define ' + var + '=' + expr)

                if self.passive: continue
                self.substitute(fileID=self.ds.fileID)

            elif self.lang.is_rgb(cmd):

                self.substitute(color=self.color)
                self.color = self.lang.get_color(cmd)
                if self.color > 30: self.cmds.append(cmd)

            else:

                self.cmds.append(cmd)

        self.substitute(color=self.color)

    def execute(self, commands=None):

        self.is_regional = True
        self.cmds  = self.get_layer_stack('cmds')

        if self.passive: return

        if commands is None:
            commands = self.cmds

        for cmd in commands:

            self.ds.reset_region()

            if self.lang.is_display(cmd):

                self.georeference()
                self.ds.pad_region(pad=8)

                expr       = []
                expression = ''.join(cmd.split(' ')[1:])

                for fld in expression.split(';'):
                    try:
                        self.fields.append(self.ds.exp(fld))
                    except:
                        self.fields.append(fld)

            elif self.lang.is_define(cmd):

                self.ds.pad_region(pad=8)
                self.ds(cmd)

            elif self.lang.is_data_service(cmd):

                self.ds(cmd)

    def get_layer_stack(self, name, layer=None):

        self.stack[name] = self.stack.get(name, {})
        stack = self.stack[name]

        if layer is not None:
            stack[layer] = stack.get(layer, [])
            return stack[layer]

        content = []
        layers  = sorted(stack.keys())

        for layer in layers:
            if layer < 0: continue
            content += stack[layer]

        for layer in layers:
            if layer != -1: continue
            content += stack[layer]

        for layer in layers:
            if layer >= -1: continue
            content += stack[layer]

        return content
            
    def georeference(self):

        qh = self.ds.query("dims")

        if qh.nt  > 1: self.is_regional = False
        if qh.nz  > 1: self.is_regional = False
        if qh.nx <= 1: self.is_regional = False
        if qh.ny <= 1: self.is_regional = False

#       Check for fixed Lat or Lon dimension

#       cmd = re.sub('`', '', cmd)
#       if not self.lang.is_latlon(cmd): return

#       values = cmd.split()[2:]
#       if len(values) <= 1:
#           self.is_regional = False
#       elif float(values[0]) == float(values[1]):
#           self.is_regional = False

        return

    def rsubstitute(self, s, **defs):

        s_interp = Template(s).safe_substitute(defs)

        if s_interp != s:
            s_interp = self.rsubstitute(s_interp, **defs)

        return s_interp

    def substitute(self, **kwargs):

        for key in kwargs:

            if key in self.special:

                val = str(kwargs[key])
                var = self.special[key]

                if val is None:
                    next

                for i in range(0,len(self.cmds)):
                    self.cmds[i] = self.cmds[i].replace(var,val)
                    self.cmds[i] = self.cmds[i].replace('&b',' ')

    def filter(self, commands, request):

        plot_only = request.get('plot_only',False)
#       no_title  = request.get('no_title',False) or plot_only
        no_label  = request.get('no_label',False) or plot_only
        no_cbar   = request.get('no_cbar',False) or plot_only

        cmds = []

        for cmd in commands:

            if cmd[0] == '#': continue

            if cmd[0] == '>':
                cmds.append(cmd[1:])
                continue

            if self.lang.is_auto(cmd): continue

            if self.lang.is_special(cmd):
                cmds.append(cmd)
                continue

            if self.lang.is_annotate(cmd) and plot_only and 'shp' not in cmd:
                continue

            if 'cbar' in cmd and no_cbar: continue

            if self.lang.is_display(cmd) and no_label:
                cmds.append('set xlab off')
                cmds.append('set ylab off')

            cmds.append(cmd)

        return cmds

    def query(self, expr):

        expr = self.eval.evaluate(field.Field(expr))
        return (self.ds.handle(), expr)

    def get_name(self):

        inum = "%03d" % (self.counter,)
        self.counter += 1
        return chr(ord('a') + self.panel) + inum

    def demand(self, **kwargs):

        self.request.update(kwargs)

    def reset(self):

        self.handle  = PlotHandle(**self.request)
        self.defs    = self.handle.__dict__

    def refine_clevs(self, clevs, nsub, type):

        clevs = [ float(clev) for clev in clevs.split() if clev != ' ' ]

        clevs_f = []
        if type.upper() == 'LOG':

            for index, clev in enumerate(clevs[0:-1]):

                dclev = ( math.log(clevs[index+1]) - math.log(clev) ) / nsub
                clevs_f.append(str(clev))

                for ns in range(2,nsub+1):
                    clev_f = math.log(clev) + (ns-1) * dclev
                    clev_f = math.exp(clev_f)
                    clev_f = "%3.2f"%clev_f
                    clevs_f.append(clev_f)

            clevs_f.append(str(clevs[-1]))

        else:

            for index, clev in enumerate(clevs[0:-1]):

                dclev = ( clevs[index+1] - clev ) / nsub
                clevs_f.append(str(clev))

                for ns in range(2,nsub+1):
                    clev_f = clev + (ns-1) * dclev
                    clevs_f.append(str(clev_f))

            clevs_f.append(str(clevs[-1]))

        return ' '.join(clevs_f)

    def refine_rgb(self, rgba, nsub):

        colors = []
        rgba_f = {}

        for color in rgba:
            colors.append([ float(col) for col in color.split() if col != ' ' ])

        nch = len(colors[0])

        for ch in range(0, nch):

            value = []
            for col in range(0, len(colors)-1):

                dch = float(colors[col+1][ch] - colors[col][ch]) / nsub
                value.append(colors[col][ch])

                if ch == 0 and nch > 4:
                    for ns in range(1, nsub):
                        value.append(value[-1] + dch)
                else:
                    for ns in range(1, nsub):
#                       value.append(int(value[-1] + dch + 0.5))
                        value.append(value[-1] + dch)

            value.append(colors[-1][ch])
            rgba_f[ch] = value

        rgba_f[nch] = [str(ch) for ch in rgba_f[0]]

        for ch in range(1, nch):
            for index, val in enumerate(rgba_f[ch]):
                rgba_f[nch][index] += ' ' + str(val)

        return rgba_f[nch]

    def refine_rgb_color(self, colors, ncols, **kwargs):

        if isinstance(colors, list):
            colors = self.to_cmap(colors)
            return self.refine_rgb_dict(colors, ncols, **kwargs)
        elif isinstance(colors, dict):
            return self.refine_rgb_dict(colors, ncols, **kwargs)
        else:
            return self.refine_rgb_name(colors, ncols)

    def to_cmap(self, clist):

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

        return cmap

    def refine_rgb_list(self, clist, ncols):

        colors = []

        for color in clist:
            rgba = [float(s)/255.0 for s in color.split() if s != ' ']
            if len(rgba) < 4: rgba.append(1.0)
            colors.append(rgba)

        cmap = mc.LinearSegmentedColormap.from_list('my_map', colors, N=ncols)

        return self.to_rgba(cmap, ncols)

    def refine_rgb_dict(self, cdict, ncols, **kwargs):

        reverse = cdict.get('reverse', 0)
        reverse = kwargs.get('reverse', reverse)
        scale   = cdict.get('scale', None)
        scale   = kwargs.get('scale', scale)
        if scale: scale = getattr(self.service, scale)

        cmap = self.unpack_cmap(cdict, **kwargs)
        cmap = SegmentedColormap('my_map', cmap, scale=scale, reverse=reverse)
        return self.to_rgba(cmap, ncols)

    def refine_rgb_name(self, cname, ncols):

        cmap = plt.get_cmap(cname)
        return self.to_rgba(cmap, ncols)

    def unpack_cmap(self, cdict, **kwargs):

        cmap = {}

        for channel in ['red', 'green', 'blue', 'alpha']:

            in_list  = cdict.get(channel, None)
            in_list  = kwargs.get(channel, in_list)

            if in_list:
                out_list = []

                for item in in_list:
                    items = item.split()
                    out_list.append([ float(c) for c in items if c != ' ' ])

                cmap[channel] = out_list

        if 'cmap' in cdict:
          cmap['cmap'] = plt.get_cmap(cdict['cmap'])

        return cmap

    def to_rgba(self, cmap, ncols):

        colors = []
        x = np.linspace(0.0, 1.0, ncols)

        for val in x:
            rgba = cmap(val)
            rgba = [ str(min(float(c)*255.0,255)) for c in rgba ]
            rgba = ' '.join(rgba)
            colors.append(rgba)

        return colors

    def set_rgb(self, rgba, nlevs, zorder=0):

        ccols = []
        for index, color in enumerate(rgba):

            c = self.lang.get_color('set rgb $* ' + color)
            ccols.append(str(c))
            self.cmd('set rgb ' + str(c) + ' ' + color, zorder=zorder)
            last_color = str(c)

        npad = nlevs - len(rgba) + 1
        cpad = '' if npad <= 0 else (last_color+' ') * npad

        self.cmd(('set ccols '+' '.join(ccols) + ' ' + cpad).strip(), zorder=zorder)

    def set_clevs(self, clevs, cmin, cmax, cint):

        if clevs: return clevs
        if cmin is None: return 
        if cmax is None: return
        if cint is None: return

        vmin = float(cmin)
        vmax = float(cmax)
        vint = float(cint)

        v    = vmin
        cout = []

        assert vmin <= vmax and vint > 0, 'invalid range for cmin/cmax/cint'
    
        while v <= vmax:
            cout.append(str(v))
            v += vint

        return ' '.join(cout)

    def set_shade(self, clevs, rgba, nsub=1, type='linear', zorder=0, **kwargs):

        if not clevs: return

        if kwargs.get('inverse', None):
            clevs = [-float(c) for c in clevs.split() if c != ' ']
            clevs = ' '.join([str(c) for c in reversed(clevs)])

        if len(clevs.split()) == 1:
            self.cmd('set rgb $* ' + rgba[0], zorder=zorder)
            self.cmd('set ccolor $*', zorder=zorder)
            return

#       rgba_f  = self.refine_rgb(rgba, nsub)
        clevs_f = self.refine_clevs(clevs, nsub, type)
        nlevs   = len(clevs_f.split())
        self.cmd('set clevs ' + clevs_f, zorder=zorder)

        rgba_f  = self.refine_rgb_color(rgba, nlevs+1, **kwargs)
        self.set_rgb(rgba_f, nlevs, zorder=zorder)

    def plot_shapes(self):

        tk    = self.config.get('toolkit', toolkit.Toolkit())
        fpath = [self.theme,'plot',self.request['field'],'shape']
        shape =  self.config(fpath, 'on')

        if str(shape) == 'off': return
        if not self.request.get('shape', 1): return

        if isinstance(shape, dict):
            map = copy.deepcopy(shape)
        else:
            map = dict(self.get_map(self.request)[0])

        map.update(self.config(['shape','default'],{}))

        for shape in map:
        
            sh_class = self.config(['shape',shape,'class'],shape)
            draw     = getattr(tk, sh_class, None)

            if draw is not None:
                attr = dict(map)
                attr.update(self.config(['shape',shape],{}))
                draw(self,map[shape],**attr)

    def plot_logos(self):

        handle  = self.handle
        request = self.request
        theme   = self.theme
        region  = request['region']
        field   = request['field']
        stream  = request['stream']

        path      = ['region', region, 'map_only']
        map_only  = self.config(path, 0)
        if map_only: return

        plot_only  = request.get('plot_only',False)
        no_logo    = request.get('no_logo',False) or plot_only
        lights_off = request.get('lights_off',False)

        if no_logo: return

        logo_path = self.config(['logo_path'],os.getcwd())

        path    = [theme, 'plot', field]
        logos   = self.config(['logo', 'default', 'layers'], [])
        logos   = self.config(['stream', stream, 'logos'], logos)
        logos   = self.config(path + ['logos'], logos)
        default = self.config(['logo', 'default'], {})
        default = { k:v for k,v in default.iteritems()
                           if isinstance(v, basestring) }

        for l in logos:

            nodes    = l.split('-')
            position = None

            if len(nodes) > 2:
                position = nodes[-1]
                name     = '-'.join(nodes[0:-1])
            else:
                name = l

            if lights_off: name += '-dark'

            logo = dict(default)
            logo.update(self.config(['logo', name], {}))
            logo.update(self.config(path + [l], {}))

            logo_file = logo['file']
            if position: logo['position'] = position

            if not os.path.isabs(logo_file):
                logo['file'] = os.path.join(logo_path, logo_file)

            self.cmd("""
              set logo file $file
              set logo position $position
              set logo size $size
              set logo margin $margin
              draw logo
              """, **logo
            )

    def plot_labels(self):

        handle   = self.handle
        request  = self.request
        region   = self.request['region']

        path      = ['region', region, 'map_only']
        map_only  = self.config(path, 0)
        if map_only: return

        plot_only = request.get('plot_only',False)
        no_title  = request.get('no_title',False) or plot_only
        no_title  = no_title or request.get('title','')

        annotate = self.config.get('annotate', {})
        path = [self.theme,'plot',request['field']]

        for ltype in self.lang.label_types:

            d = annotate.get(ltype, None)
            if not isinstance(d, dict): continue

            d_alt = self.config(path+[ltype], {})
            if isinstance(d_alt, dict):
                d.update(d_alt)
            else:
                d['string'] = d_alt

            skip_label = no_title
            if ltype in self.request: skip_label = False

            if skip_label: continue

            handle.type   = ltype
            handle.string = d.get('string',   '')
            handle.string = self.request.get(ltype, handle.string)
            handle.color  = d.get('color',    '0 0 0')
            handle.size   = d.get('size',     '0.2')
            handle.pos    = d.get('position', 'c')
            handle.margin = d.get('margin', '0.05')

            font = d.get('font', 'regular')
            if font == 'variable': font = 'regular'
            if font[0] != '$': font = '$' + font
            handle.font = font

            if not handle.string: continue

            self.cmd("""
              @TYPE $type
              @MARGIN $margin
              set font $font
              set rgb $* $color
              set string $* $pos 2 0
              set strsiz $size
              draw label $string
              """
            )

#       Draw Y-axis label if requested

        handle.ylab = self.config(path+['ylab'],   '--auto')
        self.cmd('draw ylab $ylab', zorder=-1)

    def plot_map(self):

        maps     = self.get_map(self.request)
#       basemaps = self.get_map(self.request, layers=['basemap'])
        masks    = maps[0].get('masks', [])
        masks    = self.request.get('mask', masks)
        masks    = self.get_map(self.request, layers=masks)

        self.cmd('&BASEMAP')

        self.set_map(maps)
        self.set_map_types(maps)
        self.plot_map_data(maps)
        self.plot_map_shapes(maps)
        self.plot_map_data(masks)
        self.plot_map_base(maps)

        self.cmd('&END')

#   def plot_map_labels(self):

#       map = self.get_map(self.request)[0]
#       mproj = map.get('mproj', None)
#       if mproj != 'nps' and mproj != 'sps': return

#           self.cmd("""
#             run polargrid.gs $mpvals
#             """
#           )

    def get_map(self, request, map=None, layers=None):

        """Retrieves parameters for each map layer.

        **Args:**
            request : dict : input
                Request values. Map definitions can be specified by region
                and/or field:::

                request['region']: region short name
                request['field']: field short name

            map : list : input : hidden
                List of accumulated map layers. Each map layer is a
                dictionary of parameters. This method recursively accumulates
                the map layers. The initial value should be unspecified. 

        **Returns:**
            maps: dict : output
                List of map layers. Each map layer is a dictionary of
                parameters.

        **Raises:**
            none

        **Notes:**
            none
        """

        # Overlay the map definitions in the following order: default, field
        # definition, region definition

        if map is None:

            addlayers = []
            field     = request['field']
            region    = request['region']
            map       = copy.deepcopy(self.config(['map','default'],{}))

            path = [self.theme,'plot',field,'map']
            map.update(copy.deepcopy(self.config(path,{})))
            addlayers += self.config(path+['addlayers'],[])

            path = ['region',region]
            map.update(copy.deepcopy(self.config(path)))
            addlayers += self.config(path+['addlayers'],[])

            map['mpvals'] = map.get('mpvals',None)
            map['frame']  = map.get('frame', 'on')
            map['grid']   = map.get('grid', '--auto')
            map['xlint'] = map.get('xlint', '--auto')
            map['ylint'] = map.get('ylint', '--auto')
            if not map['mpvals']: map['mpvals'] = map['lon'] + ' ' + map['lat']

            map['layers'] += addlayers

        # Recursively obtain map layers.

        maps = [map]
        if layers is None: layers = map.get('layers', [])

        for name in layers:
            layer = { k:v for k,v in map.iteritems() if k != 'layers' }
            layer.update(self.config(['map',name]))
            maps += self.get_map(request, layer)

        return maps

    def set_map(self, maps):

        self.cmd("""
          set lon $lon
          set lat $lat
          set frame $frame
          set grid $grid
          set mpdset $mpdset
          set mproj $mproj
          set mpvals $mpvals
          set xlint $xlint
          set ylint $ylint
          """, **maps[0]
        )

    def set_map_types(self, maps):

        if not maps:
            return

        self.cmd("""
          set mpt 0 off
          set mpt 1 off
          set mpt 2 off
          """
        )

        for map in maps:
            if 'mtype' in map:
                self.cmd('set rgb $* $line_color', **map)
                self.cmd('set mpt $mtype $* $line_style $line_width',**map)

    def plot_map_shapes(self, maps):

        if not maps:
            return

        shape_path = self.config(['shape_path'],os.getcwd())

        for map in maps:
            if 'shape_file' in map:
                shape_file = map['shape_file']
                if not os.path.isabs(shape_file):
                    map['shape_file'] = os.path.join(shape_path, shape_file)

                if map.get('fill_color',None):

                    self.cmd("""
                      set rgb $* $fill_color
                      set shpopts $*
                      set line $*
                      draw shp $shape_file
                      """, layer=0, **map)

                if map.get('line_color',None):

                    self.cmd("""
                      set rgb $* $line_color
                      set line $* $line_style $line_width
                      set shpopts -1
                      draw shp $shape_file
                      """, layer=1, **map)

    def plot_map_data(self, maps, zorder=None):

        if not maps: return

        handle  = self.handle
        request = self.request
        time    = request['time_dt']

        for map in maps:

            if 'expr' not in map: continue

            handle.expr    = map['expr']
            handle.gxout   = map.get('gxout', 'grfill')
            handle.gtime   = time.strftime("%Hz%d%b%Y")
            handle.name    = self.get_name()
            handle.csmooth = map.get('csmooth','off')
            handle.mask    = map.get('mask','--auto')

            if not zorder: zorder = map.get('zorder',0)
            map['zorder'] = zorder

            clevs   = map.get('clevs',None)
            cmin    = map.get('cmin',None)
            cmax    = map.get('cmax',None)
            cint    = map.get('cint',None)
            nsub    = map.get('nsub',1)
            type    = map.get('type','linear')
            cbar    = map.get('cbar','None')
            rgba    = self.config(['attribute','colorbar',cbar], cbar)
            clevs   = self.set_clevs(clevs, cmin, cmax, cint)

            if cbar in map:
                rgba = [map[cbar] for i in range(len(clevs.split()) + 1)]

            self.set_shade(clevs, rgba, nsub=nsub, type=type, zorder=zorder)

            self.cmd("""
              set dfile $#
              set time $gtime
              set lev $level
              set datawarn on
              set gxout $gxout
              set csmooth $csmooth
              set grads off
              define $name = $expr
              define $name = maskout($name,$name-$mask)
              d $name
            """, **map
            )

    def plot_map_base(self, maps, zorder=None):

        if not maps: return

        for map in maps:

            if 'service' not in map: continue

            args = json.dumps(map)
            self.cmd('>draw basemap ' + args)
              
    def get_layer(self, name):

        if self.request.get(name, 'on') == 'off': return {}

        region  = self.request['region']
        field   = self.request['field']

        path     = ['region', region, 'map_only']
        map_only = self.config(path, 0)
        if map_only: return {}

        rpath   = ['region', region, name]
        fpath   = [self.theme, 'plot', field, name]
        lpath   = [self.theme, 'layer', name]
        paths   = [lpath, fpath, rpath]

        layer = {}

#       Locate the layer in the specified hierarchy of
#       locations (region, plot, plot layer). A scalar
#       layer value with any other value than 'on' will
#       turn off the layer.

        for path in paths:
            result = self.config(path,{})
            if isinstance(result,dict):
                layer.update(result)
            elif str(result) == 'on':
                continue
            else:
                return {}

#       Resolve contour dictionaries within the layer. Each instance
#       aggregates values for one or more dependencies (i.e. by region,
#       level, time of year etc.).

        while 'cdict' in layer:
            value = self.get_attr(layer, 'cdict')
            layer = {k:v for k,v in layer.iteritems() if k != 'cdict'}
            if isinstance(value,dict): layer.update(value)

        return layer

    def get_attr(self, layer, name, default=novalue):

        region = self.request['region']
        level  = self.request['level']
        apath  = [self.theme,'attribute']

#       Return if the requested attribute name is not found.

        try:
            value = layer[name]
        except KeyError:
            if default is novalue:
                raise
            else:
                return default

#       Simply return the attribute value if it doesn't
#       reference another dependency.

#       if isinstance(value, dict): return json.dumps(value)
#       if isinstance(value, dict): return value
        if not isinstance(value, collections.Hashable): return value

        attribute = self.config(apath + [value], None)
        if attribute is None: return value

        if isinstance(attribute, list):
            return self.get_attr_dict(attribute)

        if not isinstance(attribute, dict):
            return attribute

#       Check to see if the attribute value references another
#       dependency. If so, get the value of the dependency.

        result = self.config(apath + [value, level], None)
        if result is not None: return result

        result = self.config(apath + [value, region], None)
        if result is not None: return result

        result = self.get_attr_special(attribute, None)
        if result is not None: return result

        return self.config(apath + [value, 'default'])

    def get_attr_special(self, attributes, default=None):

        time_dt = self.request['time_dt']

        for key in attributes:

            keylist = key.split('%')
            if len(keylist) <= 1: continue

            values  = keylist[0].split(',')
            tokens  = '%' + '%'.join(keylist[1:])
            vtoken  = time_dt.strftime(tokens)

            for v in values:
                if v == vtoken:
                    return attributes[key]

        return default

    def get_attr_dict(self, dlist):

        time_dt = self.request['time_dt']

        attributes = {}

        for d in dlist:

            match = True

            for key,value in d.iteritems():

                if key[0] == '%':

                    valid = time_dt.strftime(key)
                    for v in value.split(','):
                        if v == valid: break
                    else:
                        match = False

                elif key[0] == '$':

                    valid = Template(key).safe_substitute(self.request)
                    for v in str(value).split(','):
                        if v == valid: break
                    else:
                        match = False

                else:

                    continue

                if not match: break

            if match: attributes.update(d)

        return attributes

    def get_vars(self, layer):

        define    = self.get_attr(layer, 'define', '')
        var_names = [ v for v in define.split() if v != ' ' ]

        vars = {}
        for name in var_names:
            vars[name] = self.get_attr(layer, name)

        return vars

    def get_skip(self, expr):

        if self.passive: return (1, 1)

        region = self.request['region']

        fld  = self.eval.evaluate(field.Field(expr))
        fh   = self.ds.handle()
        map  = self.get_map(self.request)[0]
        lats = [ float(lat) for lat in map['lat'].split() if lat != ' ' ]
        lons = [ float(lon) for lon in map['lon'].split() if lon != ' ' ]
        mpvals = [ float(v) for v in map['mpvals'].split() if v != ' ' ]

        dx = fh.ctlinfo.dx
        dy = fh.ctlinfo.dy

        x_degrees = mpvals[1] - mpvals[0]
        y_degrees = mpvals[3] - mpvals[2]

        self.ds('set dfile ' + str(self.ds.fileID))
        self.ds('set t 1')
        self.ds('set z 1')
        self.ds('set lon ' + map['lon'])
        self.ds('set lat ' + map['lat'])
        self.ds('set mpdset ' + map['mpdset'])
        self.ds('set mproj ' + map['mproj'])
        self.ds('set mpvals ' + map['mpvals'])
        self.ds('d const(' + fld + ',0.0,-a)')
        self.ds('query gxinfo')

        x_inches  = float(self.ds.rword(3,6)) - float(self.ds.rword(3,4))
        y_inches  = float(self.ds.rword(4,6)) - float(self.ds.rword(4,4))

        x_skip = round((x_degrees / dx) / (x_inches / 0.30))
        y_skip = round((y_degrees / dy) / (y_inches / 0.30))

        i_skip = max(int(x_skip),1)
        j_skip = max(int(y_skip),1)

        return (i_skip, j_skip)

    def __iter__(self):

        index   = 0
        data    = []
        cmds    = []
        default = copy.deepcopy(self.lang.default)

        self.lang.register({'BASEMAP': not self.is_regional})

        for cmd in self.cmds:

            cmd.strip('`')
            cmds.append(cmd)
            self.lang.eval(cmd)

            if self.lang.is_display(cmd):
                n    = len(cmd.split(';'))
                data = self.fields[index:index+n]
                index += n

            if self.lang.is_annotate(cmd):
                data = self.lang.value

            if self.lang.is_action(cmd):
                state   = self.lang.state
                macro   = self.lang.macro
                yield PlotObject(state, macro, data, cmds, default)

                cmds = []

class PlotHandle(object):

    def __init__(self, **kwargs):

        self.__dict__.update(kwargs)

class PlotObject(object):

    def __init__(self, state, macro, data=None, cmds=None, default=None):

        self.state   = state
        self.macro   = macro
        self.data    = data
        self.cmds    = cmds
        self.default = default
