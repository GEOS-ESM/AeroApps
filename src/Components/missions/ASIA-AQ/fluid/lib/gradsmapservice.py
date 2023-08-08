import sys
import uuid
import json
import math
import re
import os
import io

from string import *

import numpy.ma as ma
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
from PIL import Image, ImageChops, ImageEnhance, ImageDraw, ImageFont
from PIL.ImageColor import getcolor, getrgb
from PIL.ImageOps import grayscale

from mapservice import *

novalue = object()

class Service(MapService):

    def __init__(self, *args, **kwargs):

        super(Service,self).__init__(*args, **kwargs)

        self.map    = None
        self.cbar   = None
        self.logos  = []
        self.labels = {}

        self.handler = {'DRAW_LINE'       : self.draw_line,
                        'DRAW_MARK'       : self.draw_mark,
                        'DRAW_POLYF'      : self.draw_polyf,
                        'DRAW_REC'        : self.draw_rec,
                        'DRAW_RECF'       : self.draw_recf,
                        'DRAW_WXSYM'      : self.default,
                        'DRAW_HILO'       : self.draw_hilo,
                        'DRAW_SHP'        : self.default,
                        'DRAW_STRING'     : self.draw_string,
                        'DRAW_SUBTITLE'   : self.draw_subtitle,
                        'DRAW_TMSTRING'   : self.draw_tmstring,
                        'RUN'             : self.default,
                        'DRAW_LOGO'       : self.save_logo,
                        'DRAW_LABEL'      : self.save_label,
                        'DRAW_CBAR'       : self.draw_cbar,
                        'DRAW_ARROW'      : self.draw_arrow,
                        'DRAW_GRID'       : self.draw_grid,
                        'DRAW_BASEMAP'    : self.draw_basemap,
                        'BASEMAP'         : self.ignore,
                        'DISPLAY_GRID'    : self.display_grid,
                        'DISPLAY_VECTOR'  : self.display_vector,
                        'DISPLAY_SHADED'  : self.display,
                        'DISPLAY_CONTOUR' : self.display,
                        'DISPLAY_GRFILL'  : self.display,
                        'DISPLAY_BARB'    : self.display,
                        'DISPLAY_STREAM'  : self.display
                       }

    def get_capabilities(self, request):
        return None

    def get_maps(self, plots):

        if len(plots) == 1: return self.get_map(plots[0])

        for i, plot in enumerate(plots):
            path = os.path.dirname(plot.request['oname'])
            name = str(os.getpid()) + '-%d.png'%(i,)
            pathname = os.path.join(path,name)
            self.get_map(plot, oname=pathname)

        return self.aggregate(plots)

    def aggregate(self, plots):

        cbar    = ()
        hotspot = []
        plot    = plots[0]
        frames  = self.get_frames(plots)

        geometry = plot.options.get('geometry','1024x768').split('x')
        geometry = [int(s) for s in geometry]

        tight = plot.request.get('tight', False)

        if 'bbox' in plot.request:
            bbox = plot.request['bbox']
            bbox = [int(x) for x in bbox.split()]

        color = (255, 255, 255)
        if plot.options.get('lights_off', False): color = (0, 0, 0)

        bg = Image.new('RGB', geometry, color)

        for i, plot in enumerate(plots):

            xfc, yfc, pos = frames[i]

            gp = plot.request.get('geometry','1024x768')
            offset = re.match(r'.+([+-]\d+)([+-]\d+)',gp)

            xshift, yshift = (0, 0)

            if offset:
                xshift = int(offset.group(1))
                yshift = int(offset.group(2))

            path = os.path.dirname(plot.request['oname'])
            name = str(os.getpid()) + '-%d.png'%(i,)
            file = os.path.join(path,name)

            im   = self.image_trim(Image.open(file))
            xc   = int(im.size[0] / 2.0)
            yc   = int(im.size[1] / 2.0)

            if pos == 'cbar':
                cbar = (file, im, xfc-xc+xshift, yfc-yc+yshift-im.size[1])
            else:
                bg.paste(im, (xfc-xc+xshift,yfc-yc+yshift))
                im.close()
                os.remove(file)
                url = plot.request.get('url', None)
                if url: hotspot.append({'x':xfc+xshift, 'y':yfc+yshift, 'r':xc, 'url':url})

      # Draw text on image.

        W  = geometry[0]
        H  = geometry[1]
        im = Image.new(bg.mode, bg.size, bg.getpixel((0,0)))

        if 'bbox' not in plot.request or tight:
            diff = ImageChops.difference(bg, im)
            diff = ImageChops.add(diff, diff, 2.0, -100)
            bbox = diff.getbbox()
            im.close()

        color = (0, 0, 0)
        if plot.options.get('lights_off', False): color = (255, 255, 255)

        d = ImageDraw.Draw(bg)

        f = ImageFont.truetype(self.config['bfont'], 36)
        s = self.get_label('title')
        w, h1 = d.textsize(s, font=f)
        d.text( ((W-w)/2,bbox[1]-h1-5), s, font=f, fill=color)

        f = ImageFont.truetype(self.config['bfont'], 16)
        s = self.get_label('header')
        w, h2 = d.textsize(s, font=f)
        d.text( ((W-w)/2,bbox[1]-h1-h2-5), s, font=f, fill=color)

#       f = ImageFont.truetype(self.config['bfont'], 36)
#       s = self.get_label('main')
#       w, h2 = d.textsize(s, font=f)
#       d.text( ((W-w)/2,(h1+h2)/2+5), s, font=f, fill=color)

        f = ImageFont.truetype(self.config['bfont'], 24)
        s = self.get_label('xlabel')
        w, h = d.textsize(s, font=f)
        d.text( ((W-w)/2,bbox[3]+5), s, font=f, fill=color)

      # Draw the colorbar if it exists

        if cbar: 

            im = Image.new(bg.mode, bg.size, bg.getpixel((0,0)))
            diff = ImageChops.difference(bg, im)
            diff = ImageChops.add(diff, diff, 2.0, -100)
            bbox = diff.getbbox()
            im.close()

            file, im, xc, yc = cbar
            xc = (W - im.size[0]) / 2
            yc = bbox[3] + 10
            bg.paste(im, (xc,yc))
            im.close()
            os.remove(file)

      # Save image

        oname = plots[0].request['oname']
        bg.save(oname)

        bg.close()

      # Write out navigational information

        if self.request.get('navigate', 'on') != 'off':

            navfile = '.'.join(oname.split('.')[0:-1]) + '.nav'
            if hotspot:
                with open(navfile, 'w') as f: json.dump(hotspot, f)

        return oname

    def get_label(self, type):

        request = self.options

        annotate = self.config.get('annotate', {})
        annotate = annotate.get(type, {})

        path  = [self.theme,'plot',request['field']]
        label = annotate.get('string', '')
        label = self.config(path+[type], label)

        if isinstance(label, dict):
            str = label.get('string', '')
        else:
            str = label

        str = Template(str).safe_substitute(self.defs)
        str = Template(str).safe_substitute(self.defs)

        return str

    def image_trim(self, im):

        bg = Image.new(im.mode, im.size, im.getpixel((0,0)))
        diff = ImageChops.difference(im, bg)
        diff = ImageChops.add(diff, diff, 2.0, -100)
        bbox = diff.getbbox()
        if bbox:
            return im.crop(bbox)

    def get_map(self, plot, **kwargs):

        self.logos   = []
        self.labels  = {}
        self.cbar    = None
        self.ccols   = None
        self.clevs   = None
        self.maps    = []
        self.theme   = plot.theme
        self.request = plot.request
        self.options = plot.options
        self.config  = plot.config
        self.defs    = plot.defs
        self.slice   = False
        self.dims    = None
        self.imap    = []
        self.oname   = kwargs.get('oname', plot.request['oname'])

        for obj in plot:

            handler = self.handler.get(obj.macro, self.default)
            handler(obj)

            self.ccols = obj.state.get('ccols', self.ccols)
            self.clevs = obj.state.get('clevs', self.clevs)

        return self.imshow()

    def imshow(self):

        cbar_only = self.request.get('cbar_only', False)

        if cbar_only: self.clear_map()

        if not cbar_only: self.make_labels()
        if self.cbar: self.cbar.draw(**self.request)

        background = self.draw_map(zorder=0)

        geometry = self.get_geometry()
        geometry = 'x' + str(geometry[0]) + ' ' + 'y' + str(geometry[1])

        img = self.oname

        t_color = 1
        if self.request.get('lights_off', False): t_color = 0

        if background:
            self.ds('gxprint ' + img + ' ' + geometry + ' -b ' +
                        background + ' -t  ' + str(t_color))
        else: 
            self.ds('gxprint ' + img + ' ' + geometry)

        if not cbar_only: self.draw_logo(img, self.logos)

        self.navigate(img)

        return img

    def navigate(self, img):

        if not self.imap: return

        path  = [self.theme,'plot',self.request['field']]
        if self.config(path+['navigate'], 'on') == 'off': return
        if self.request.get('navigate', 'on') == 'off': return

        data    = []
        ax      = Axes(self.ds)
        navfile = '.'.join(img.split('.')[0:-1]) + '.nav'

      # Get the image size.

        im    = Image.open(img)
        xsize = im.size[0]
        ysize = im.size[1]
        im.close()

      # Extract navigational information in
      # pixel space.

        for imap in self.imap:

            xinch  = imap[0]
            yinch  = imap[1]
            radius = imap[2]
            hover  = ' '.join(imap[3:-1])
            url    = imap[-1]
            xp     = int(xinch / 11.0 * xsize)
            yp     = ysize - int(yinch /  8.5 * ysize)
            rp     = int(radius / 11.0 * xsize)

            data.append({'x':xp, 'y':yp, 'r':rp, 'hover':hover, 'url':url})

      # Write out navigational information.

        with open(navfile, 'w') as f: json.dump(data, f)

    def clear_map(self):

        ax = Axes(self.ds)

        lights_off = self.request.get('lights_off', False)

        if lights_off:
            self.ds('set rgb 30 0 0 0')
        else:
            self.ds('set rgb 30 255 255 255')

        self.ds('set line 30')
        self.ds('draw recf %f %f %f %f'%(ax.LL+ax.UR))

    def get_geometry(self):

        g = self.request.get('geometry','1024x768')
        g = re.match(r'(\d+)x(\d+)', g)
        g = (int(g.group(1)), int(g.group(2)))

        gs = self.options.get('geometry','1024x768')
        gs = re.match(r'(\d+)x(\d+)', gs)
        gs = (int(gs.group(1)), int(gs.group(2)))

#       g = self.request.get('geometry','1024x768').split('x')
#       g = [int(s) for s in g]

#       gs = self.options.get('geometry','1024x768').split('x')
#       gs = [int(s) for s in gs]

        xsize = int((float(g[0]) / 1024.0) * float(gs[0]))
        ysize = int((float(g[1]) /  768.0) * float(gs[1]))

        return (xsize, ysize)

    def get_frames(self, plots):

        gs = self.options.get('geometry','1024x768').split('x')
        gs = [int(s) for s in gs]

        bbox = plots[0].request.get('bbox','0 60 1024 728')
        bbox = [int(x) for x in bbox.split()]

        nx = 0
        ny = 0

        for plot in plots:

            pos = str(plot.request.get('pos', 11))

            if pos.isdigit():
                nx = max(int(pos[0]),  nx)
                ny = max(int(pos[-1]), ny)

        layout = plots[0].request.get('layout', None)
        if layout:
            nx = int(layout.split('x')[0])
            ny = int(layout.split('x')[-1])

        xs    = bbox[0]
        xe    = bbox[2]
        ys    = bbox[1]
        ye    = bbox[3]

        xsize = int( (xe-xs) / float(nx))
        ysize = int( (ye-ys) / float(ny))
        xr    = int(xsize / 2.0)
        yr    = int(ysize / 2.0)

        i = 0
        j = 0
        frames = []

        for plot in plots:

            pos = str(plot.request.get('pos', 11))

            if pos.isdigit():

                i = int(pos[0]) - 1
                j = int(pos[-1]) - 1
                frames.append( (i*xsize+xr+xs, j*ysize+yr+ys, pos) )
                i = (i+1)%nx
                if i == 0: j = j + 1
            else:
                frames.append( (gs[0]/2, gs[1], pos) )
#               frames.append( ((bbox[2]-bbox[0])/2, bbox[3], pos) )

        return frames

    def draw_grid(self, obj):

        return

#       Query key parameters

        lonmin = float(obj.state['lon'].split()[0])
        lonmax = float(obj.state['lon'].split()[1])
        latmin = float(obj.state['lat'].split()[0])
        latmax = float(obj.state['lat'].split()[1])

        midlat     = (latmax+latmin)/2.0
        quarterlat = (latmax+midlat)/2.0

        self.ds('query dims')

        lonstart = float(self.ds.rword(2,6))
        lonend   = float(self.ds.rword(2,8))
        midlon   = (lonstart+lonend)/2.0

    def display_grid(self, obj):

        cmds   = [ cmd for cmd in obj.cmds if self.is_valid(cmd) ]
        cmark  = obj.state.get('cmark', '0')
        digsiz = obj.state.get('digsiz', '0.1')

        if cmark == '0':
            cmds = '\n'.join(cmds)
            self.ds(cmds)
            return

        cmds = '\n'.join(cmds[0:-1])
        self.ds(cmds)

        ax    = Axes(self.ds)
        data  = obj.data[0]
        grid  = data.grid
        kdim, jdim, idim = ma.shape(obj.data)

        for j in range(0,jdim):
            for i in range(0,idim):

                if not data[j,i]: continue

                lon = grid.lon[i]

                x, y = self.w2xy(lon, grid.lat[j])
                if not ax.is_clipped((x,y)):
                    self.ds('draw mark %s %f %f %s'%(cmark,x,y,digsiz))
                    continue

                if lon >= 0.0 and lon <= 180.0: continue

                if grid.lon[i] < 0.0: lon += 360.0
                if grid.lon[i] > 180.0: lon -= 360.0

                x, y = self.w2xy(lon, grid.lat[j])
                if not ax.is_clipped((x,y)):
                    self.ds('draw mark %s %f %f %s'%(cmark,x,y,digsiz))

    def display_slice(self, obj):

        cmd  = obj.cmds[-1].split()[1:]
        vars = ''.join(cmd).split(';')

        cmds = [ cmd for cmd in obj.cmds if self.lang.is_file(cmd) ]
        cmds = cmds + [ cmd for cmd in obj.cmds if self.lang.is_dimension(cmd) ]
        cmds = '\n'.join(cmds)
        self.ds(cmds)

        qh = self.ds.query("dims")
        if qh.z_state == 'fixed' and self.slice: return self.slice_2D(obj)

        if qh.x_state == 'fixed': return
        if qh.y_state == 'fixed': return
        if qh.z_state == 'fixed': return

        self.slice = True
        self.dims  = qh

        return self.slice_3D(obj)

    def slice_3D(self, obj):

      # Execute the commands up to but not including
      # the display command.
      # --------------------------------------------

        cmds = [ cmd for cmd in obj.cmds if self.is_valid(cmd) ]
        cmds = '\n'.join(cmds[0:-1])

        self.ds(cmds)

      # Extract the variable names from the
      # display command (more than one for vectors).
      # --------------------------------------------

        cmd  = obj.cmds[-1].split()[1:]
        vars = ''.join(cmd).split(';')

      # Initialize the dimension environment for
      # collecting profiles.
      # ----------------------------------------

        lon1, lon2 = self.dims.lon
        x1, x2     = self.dims.x
        lat1, lat2 = self.dims.lat
        y1, y2     = self.dims.y

        slice = self.get_state(obj.state,'SLICE').split()
        if slice: lon1,lat1,lon2,lat2 = tuple([float(v) for v in slice])

        self.ds('set x %d'%(int(x1),))
        self.ds('set y %d'%(int(y1),))

      # Free up the collection grids
      # (one grid per variable).
      # ----------------------------

        for index, var in enumerate(vars):
            self.ds('collect %d free'%(index+1,))

      # Collect the profiles by sampling along
      # the slice.
      # --------------------------------------

        skip     = (1,1)
        template = 'collect %d gr2stn(%s,%f,%f)'

        lon, lat = (lon1, lat1)
        while lon <= lon2:
            lat = lat1 + (lat2-lat1)*(lon-lon1) / (lon2-lon1)

            for index, var in enumerate(vars):
                match = re.match('skip\(\w+,(\d+),(\d+)\)', var)
                if match: skip = (match.group(1), match.group(2))
                var = var.split('(')[-1].split(',')[0]
                self.ds(template%(index+1, var, lon, lat))

            lon = lon + 0.1

      # Execute the display command using the 
      # collected profile grids in place of the
      # original variables.
      # ---------------------------------------

        vars    = ['coll2gr(%d,-u)'%(i,) for i in range(1,len(vars)+1)]
        tmpl    = 'skip(###,%s,%s)'%skip
        vars[0] = tmpl.replace('###',vars[0])

        self.ds('set x %d %d'%(int(x1), int(x2)))
        self.ds('d ' + ';'.join(vars))

        obj.cmds = []

    def slice_2D(self, obj):

        cmds = [ cmd for cmd in obj.cmds if self.is_valid(cmd) ]
        cmds = '\n'.join(cmds[0:-1])

        self.ds(cmds)

        var  = obj.cmds[-1].split()[-1]

        lon1, lon2 = self.dims.lon
        x1, x2     = self.dims.x
        lat1, lat2 = self.dims.lat
        y1, y2     = self.dims.y
        lev1, lev2 = self.dims.lev
        z1, z2     = self.dims.z

        slice = self.get_state(obj.state,'SLICE').split()
        if slice: lon1,lat1,lon2,lat2 = tuple([float(v) for v in slice])

        vertices = []
        lon, lat = (lon1, lat1)

        while lon <= lon2:

            lat = lat1 + (lat2-lat1)*(lon-lon1) / (lon2-lon1)

            self.ds('set gxout print')
            self.ds('set lat %f'%(lat,))
            self.ds('set lon %f'%(lon,))
            self.ds('set lev %f'%(lev1,))
            self.ds('d ' + var)

            z = self.ds.rword(2,1)
            x, z = self.w2xy(lon, z)
            vertices.append( (x, z) )


            lon = lon + 0.1

        self.ds('set lon %f %f'%(lon1, lon2))
        self.ds('set lev %f %f'%(lev1, lev2))

        p1 = vertices[0]
        for p in vertices[1:]:
            self.ds('draw line %s %s %s %s'%(p1+p))
            p1 = p

        obj.cmds = []

    def display_vector(self, obj):

        cmds = [ cmd for cmd in obj.cmds if self.is_valid(cmd) ]

        arrfill   = self.get_state(obj.state,'ARRFILL',False)
        arrscl    = self.get_state(obj.state,'arrscl').split()
        arrowhead = self.get_state(obj.state,'arrowhead').split()

        if not arrfill:
            cmds = '\n'.join(cmds)
            self.ds(cmds)
            return

        arhdwd  = 0.08
        arhdln  = 0.08
        vsize   = 0.5
        refspd  = 1.0

        if len(arrscl) >= 1: vsize  = float(arrscl[0])
        if len(arrscl) >= 2: refspd = float(arrscl[1])
        if arrowhead:
            arhdwd = float(arrowhead[0])
            arhdln = float(arrowhead[0])

        cmds = '\n'.join(cmds[0:-1])
        self.ds(cmds)

        ax    = Axes(self.ds)
        uwnd  = obj.data[0]
        vwnd  = obj.data[1]
        grid  = uwnd.grid
        kdim, jdim, idim = ma.shape(obj.data)

        for j in range(0,jdim):
            for i in range(0,idim):

                if uwnd[j,i] and vwnd[j,i]:

                    x, y = self.w2xy(grid.lon[i], grid.lat[j])
                    if x < ax.xlow or x > ax.xhigh: continue
                    if y < ax.ylow or y > ax.yhigh: continue

                    rlat   = grid.lat[j] * math.pi / 180.0

                    u1     = uwnd[j,i]
                    v1     = vwnd[j,i]
                    s1     = math.sqrt(u1**2 + v1**2)

                    u2     = u1 / math.cos(rlat)
                    v2     = v1
                    s2     = math.sqrt(u2**2 + v2**2)
                    phi2   = math.asin(v2/s2)

                    phi3   = phi2
                    s3     = s1
                    u3     = s3 * math.cos(phi3)
                    v3     = s3 * math.sin(phi3)

                    uadj   = math.copysign(u3,u1)
                    vadj   = math.copysign(v3,v1)

                    uadj = u1
                    vadj = v1

#                   xb     = x - (vsize * uadj/refspd)
#                   xa     = x + (vsize * uadj/refspd)
#                   yb     = y - (vsize * vadj/refspd)
#                   ya     = y + (vsize * vadj/refspd)

                    xb     = x
                    xa     = x + (vsize * uadj/refspd)
                    yb     = y
                    ya     = y + (vsize * vadj/refspd)
                    vlen   = math.sqrt((xa-xb)**2+(ya-yb)**2)

                    xc     = xa - (arhdln / vlen) * (xa - xb)
                    yc     = ya - (arhdln / vlen) * (ya - yb)

                    xd     = xc - (0.5 * arhdwd/arhdln) * (ya - yc)
                    xe     = xc + (0.5 * arhdwd/arhdln) * (ya - yc)
                    yd     = yc + (0.5 * arhdwd/arhdln) * (xa - xc)
                    ye     = yc - (0.5 * arhdwd/arhdln) * (xa - xc)

                    if xa < ax.xlow or xa > ax.xhigh: continue
                    if xb < ax.xlow or xb > ax.xhigh: continue
                    if ya < ax.ylow or ya > ax.yhigh: continue
                    if yb < ax.ylow or yb > ax.yhigh: continue

                    self.ds('draw line %f %f %f %f'%(xa,ya,xb,yb))
                    self.ds('draw line %f %f %f %f'%(xa,ya,xd,yd))
                    self.ds('draw line %f %f %f %f'%(xa,ya,xe,ye))
                    self.ds('draw line %f %f %f %f'%(xd,yd,xe,ye))
                    self.ds('draw polyf %f %f %f %f %f %f'%(xa,ya,xd,yd,xe,ye))

    def draw_hilo(self, obj):

        argstr = obj.cmds[-1][10:]
        args   = json.loads(argstr)

        self.make_hilo('l', **args)
        self.make_hilo('h', **args)

    def make_hilo(self, sense, **kwargs):

        ax     = Axes(self.ds)

        expr   = kwargs['expr']
        method = kwargs.get('algo', 'CL')
        radius = str(kwargs.get('radius', 1000))
        cint   = str(kwargs.get('cint', 300))
        color  = str(kwargs.get('color', '0 0 0'))

        index  = self.lang.get_color('set rgb $* ' + color)

        cmd = ' '.join(['mfhilo',expr,method,sense,radius,cint])
#       cmd = ' '.join(['mfhilo',expr,'GR',sense,'d 30 16.8 129.7'])
        self.ds(cmd)

        sense = sense.upper()
        mark = kwargs['hmark']
        if sense == 'L': mark = kwargs['lmark']

        i = 2
        points = []

        while self.ds.rword(i,1) == sense:

            lat = self.ds.rword(i,2)
            lon = self.ds.rword(i,3)
            val = self.ds.rword(i,5)
            if ('e' not in val): points.append( (lon,lat,val) )

            i = i + 1

        for (lon, lat, val) in points:

            x, y = self.w2xy(lon, lat)
            if x < ax.xlow+0.2 or x > ax.xhigh-0.2: continue
            if y < ax.ylow+0.2 or y > ax.yhigh-0.2: continue

            self.ds('set rgb %d %s'%(index, color))
            self.ds('set strsiz 0.28')
            self.ds('set string %d c 9'%(index,))
            self.ds('draw string %f %f %s'%(x, y, mark))
            self.ds('set strsiz 0.15')
            self.ds('set string %d c 9'%(index,))
            self.ds('draw string %f %f %s'%(x, y-0.25, int(round(float(val)))))

    def draw_cbar(self, obj):

        colors = [c for c in self.ccols.split() if c != ' ']
        levels = [l for l in self.clevs.split() if l != ' ']

        args    = re.sub("\s+"," ",obj.cmds[-1].strip())
        options =  json.loads(args[9:])

        if not self.cbar:

            ax = Axes(self.ds)
            ratio = (ax.xhigh - ax.xlow) / (ax.yhigh - ax.ylow)

            cbar_only = self.request.get('cbar_only', False)
            if cbar_only: ratio = 2.0

            if ratio < 1.5:
                self.cbar = VerticalColorbar(self.ds)
            else:
                self.cbar = HorizontalColorbar(self.ds)

        self.cbar.add(colors, levels, **options)

    def draw_arrow(self, obj):

        ax    = Axes(self.ds)
        args  = re.sub("\s+"," ",obj.cmds[-1].strip()).split()

        onoff  = args[2]
        if onoff != 'on': return

        veclen = float(args[3])
        scale  = args[4]
        pos    = ax.get(args[5], ax.Lr)
        
        x = pos[0] - veclen / 2.0
        y = pos[1] + 0.2

        xv1   = x   + veclen / 2.0
        xv2   = x   - veclen / 2.0
        xa    = xv1 - 0.05
        ya1   = y   + 0.025
        ya2   = y   - 0.025
        ys    = y   - 0.1

        self.ds('set line 1 1 4')
        self.ds('draw line %f %f %f %f'%(xv1,y,xv2,y))
        self.ds('draw line %f %f %f %f'%(xa,ya1,xv1,y))
        self.ds('draw line %f %f %f %f'%(xa,ya2,xv1,y))
        self.ds('set font 13')
        self.ds('set string 1 c')
        self.ds('set strsiz 0.15')
        self.ds('draw string %f %f %s'%(x,ys,scale))

    def display(self, obj):

        self.display_slice(obj)
        self.default(obj)
    
    def default(self, obj):

        cmds = [ cmd for cmd in obj.cmds if self.is_valid(cmd) ]
        cmds = '\n'.join(cmds)

        self.ds(cmds)

    def ignore(self, obj): pass

    def save_logo(self, obj):

        cmds = [ cmd.split() for cmd in obj.cmds if self.is_valid(cmd) ]
        params = { c[2]:c[3] for c in cmds if len(c) >= 4 }

        self.logos.append(params)

    def save_label(self, obj):

        cmds      = [ cmd for cmd in obj.cmds if self.is_valid(cmd) ]
        drawcmds  = [ cmd for cmd in cmds if cmd[0:4] == 'draw' ]

        if not drawcmds: return

        attr    = [ c.split() for c in obj.cmds if self.lang.is_attribute(c)]
        attr    = { c[0][1:]:c[1:] for c in attr }

        type = attr['TYPE'][0]
        self.labels[type] = obj.cmds

    def draw_logo(self, to_img, logos):

#       Get the bounding box of the background image after trimming unused
#       space on the edges.

        bg   = Image.open(to_img)
        mask = Image.new(bg.mode, bg.size, bg.getpixel((0,0)))
        diff = ImageChops.difference(bg, mask)
        bbox = diff.getbbox()
        mask.close()

#       Compute the GrADS frame in pixel space.

        self.ds('query gxinfo')

        xlo  = int(float(self.ds.rword(3,4)) / 11.0 * bg.size[0])
        xhi  = int(float(self.ds.rword(3,6)) / 11.0 * bg.size[0])
        ylo  = int(float(self.ds.rword(4,4)) / 8.5  * bg.size[1])
        yhi  = int(float(self.ds.rword(4,6)) / 8.5  * bg.size[1])
        y0   = bg.size[1] - yhi
        y1   = bg.size[1] - ylo

        frame = (xlo, y0, xhi, y1)

#       Set the bounding box to be at least 1-inch above or below
#       the GrADS frame.

        margin     = int((1.0 / 8.5)  * bg.size[1])
        top_margin = int(max(frame[1] - bbox[1] + margin/2, margin))
        bot_margin = int(max(bbox[3] - frame[3] + margin/2, margin))

        y0 = max(0, frame[1] - top_margin)
        y1 = min(bg.size[1], frame[3] + bot_margin)

        bbox = (0, y0, bg.size[0], y1)

#       Paste each logo onto the background image.

        for logo in logos:

            file    = logo.get('file', None)
            pos     = logo.get('position', 'ul')
            size    = float(logo.get('size', '10'))
            margin  = float(logo.get('margin', '1'))

#           Re-size the logo image and maintain the aspect ratio.

            img     = Image.open(file)
            ratio   = float(img.size[1]) / float(img.size[0])

            xsize   = int(size / 100.0 * bg.size[0])
            ysize   = int(xsize * ratio)
            xmargin = int(margin / 100.0 * bg.size[0])
            ymargin = int(margin / 100.0 * bg.size[1])
            xpos    = xmargin
            ypos    = ymargin

            img     = img.resize((xsize, ysize), Image.ANTIALIAS)

#           Paste the logo onto the background image.

            if pos == 'ul':

                xpos = bbox[0] + xmargin
                ypos = bbox[1] + ymargin

            elif pos == 'ur':

                xpos = bbox[2] - img.size[0] - xmargin
                ypos = bbox[1] + ymargin

            elif pos == 'll':

                xpos = bbox[0] + xmargin
                ypos = bbox[3] - img.size[1] - ymargin

            elif pos == 'lr':

                xpos = bbox[2] - img.size[0] - xmargin
                ypos = bbox[3] - img.size[1] - ymargin

            bg.paste(img, (xpos, ypos), img)
            img.close()

#       Save the background image with logos.

        bg.save(to_img, "PNG")
        bg.close()

    def draw_line(self, obj):

        clip = int(self.get_state(obj.state,'CLIP',1))
        ax   = Axes(self.ds,clip=clip)

        cmds = [ cmd for cmd in obj.cmds[0:-1] if self.is_valid(cmd) ]
        cmds = '\n'.join(cmds)
        self.ds(cmds)

        cmd  = re.sub("\s+"," ",obj.cmds[-1].strip()).split()

        x1, y1 = self.w2xy(cmd[2], cmd[3])
        x2, y2 = self.w2xy(cmd[4], cmd[5])

        line   = ax.clip([(x1,y1), (x2,y2)])
        x1, y1 = line[0]
        x2, y2 = line[1]

        self.ds('draw line %s %s %s %s'%(x1,y1,x2,y2))

    def draw_string(self, obj):

        clip = int(self.get_state(obj.state,'CLIP',1))
        ax   = Axes(self.ds,clip=clip)

        cmds = [ cmd for cmd in obj.cmds[0:-1] if self.is_valid(cmd) ]
        cmds = '\n'.join(cmds)
        self.ds(cmds)

        cmd  = re.sub("\s+"," ",obj.cmds[-1].strip()).split()

        x, y = self.w2xy(cmd[2], cmd[3])
        text = ' '.join(cmd[4:])

        if not ax.is_clipped((x,y)):
            self.ds('draw string %s %s %s'%(x,y,text))

    def draw_rec(self, obj):

        cmds = [ cmd for cmd in obj.cmds[0:-1] if self.is_valid(cmd) ]
        cmds = '\n'.join(cmds)
        self.ds(cmds)

        cmd  = re.sub("\s+"," ",obj.cmds[-1].strip()).split()

        x1, y1 = self.w2xy(cmd[2], cmd[3])
        x2, y2 = self.w2xy(cmd[4], cmd[5])

        self.ds('draw rec %s %s %s %s'%(x1,y1,x2,y2))

    def draw_recf(self, obj):

        cmds = [ cmd for cmd in obj.cmds[0:-1] if self.is_valid(cmd) ]
        cmds = '\n'.join(cmds)
        self.ds(cmds)

        cmd  = re.sub("\s+"," ",obj.cmds[-1].strip()).split()

        x1, y1 = self.w2xy(cmd[2], cmd[3])
        x2, y2 = self.w2xy(cmd[4], cmd[5])

        self.ds('draw recf %s %s %s %s'%(x1,y1,x2,y2))

    def draw_polyf(self, obj):

        cmds = [ cmd for cmd in obj.cmds[0:-1] if self.is_valid(cmd) ]
        cmds = '\n'.join(cmds)
        self.ds(cmds)

        cmd  = re.sub("\s+"," ",obj.cmds[-1].strip()).split()[2:]

        for i in range(0,len(cmd)-1,2):
            cmd[i], cmd[i+1] = self.w2xy(cmd[i], cmd[i+1])
            if i > 0: self.ds('draw line %s %s %s %s'%(cmd[i-2],cmd[i-1],cmd[i],cmd[i+1]))

#       print ' '.join(cmd)
#       self.ds('draw polyf ' + ' '.join(cmd))

    def draw_mark(self, obj):

        ax    = Axes(self.ds)

        cmds = [ cmd for cmd in obj.cmds[0:-1] if self.is_valid(cmd) ]
        cmds = '\n'.join(cmds)
        self.ds(cmds)

        cmd  = re.sub("\s+"," ",obj.cmds[-1].strip()).split()

        x, y = self.w2xy(cmd[3], cmd[4])

        if x < ax.xlow or x > ax.xhigh: return
        if y < ax.ylow or y > ax.yhigh: return

        self.ds('draw mark %s %f %f %s'%(cmd[2],x,y,cmd[5]))

        if len(cmd) > 6: self.imap.append((x,y,float(cmd[5])) + tuple(cmd[6:]))

    def draw_basemap(self, obj):

        argstr = obj.cmds[-1][13:]
        bmap   = json.loads(argstr)

        bmap['lat']    = [float(v) for v in bmap['lat'].split() if v != ' ']
        bmap['lon']    = [float(v) for v in bmap['lon'].split() if v != ' ']
        bmap['mpvals'] = [float(v) for v in bmap['mpvals'].split() if v != ' ']

#       geometry           = self.request.get('geometry','1024x768').split('x')
        bmap['geometry']   = self.get_geometry()
#       bmap['geometry']   = [int(v) for v in geometry]
        bmap['lights_off'] = self.request.get('lights_off', False)

        color_keys = ['land_color', 'water_color', 'ocean_color',
                      'lake_color', 'line_color', 'tint_color',
                      'land_tint_color', 'water_tint_color']

        for key in color_keys:
            color     = bmap.get(key, None)
            if not color: color = '0 0 0 0'
            bmap[key] = tuple([float(c)/255.0 for c in color.split() 
                                                           if c != ' '])

        if self.is_clear(bmap['land_tint_color']):
            bmap['land_tint_color'] = bmap['tint_color']

        if self.is_clear(bmap['water_tint_color']):
            bmap['water_tint_color'] = bmap['tint_color']

        default                  = float(bmap.get('brightness', '1.0'))
        bmap['brightness']       = default
        bmap['land_brightness']  = float(bmap.get('land_brightness', default))
        bmap['water_brightness'] = float(bmap.get('water_brightness', default))

        bmap['face_color'] = 'white'
        if bmap['lights_off']: bmap['face_color'] = 'black'

        bmap['grayscale'] = bmap.get('grayscale', False)

        self.maps.append(bmap)

    def draw_map(self, zorder=0):

        bmaps = [b for b in self.maps if b.get('zorder',0) == zorder]
        if not bmaps: return None

        name  = self.get_map_key(bmaps) + '.png'
        path  = self.config.get('map_path', os.getcwd())
        path  = Template(path).safe_substitute(os.environ)
        fname = os.path.join(path, name)
        if os.path.isfile(fname): return fname

        try:
            os.makedirs(path, 0755)
        except:
            pass

        m = self.get_map_handle(**bmaps[0])

        self.create_mask(m, **bmaps[0])

        for bmap in bmaps: self.add_map_layer(m, bmap)

        self.save_map(fname, **bmaps[0])

        return fname

    def save_map(self, fname, **kwargs):

        face_color        = kwargs['face_color']
        isGrayscale       = kwargs['grayscale']
        brightness        = kwargs['brightness']
        land_brightness   = kwargs['land_brightness']
        water_brightness  = kwargs['water_brightness']
        water_color       = [int(c * 255) for c in kwargs['water_color']]
        land_color        = [int(c * 255) for c in kwargs['land_color']]
        tint_color        = [int(c * 255) for c in kwargs['tint_color']]
        land_tint_color   = [int(c * 255) for c in kwargs['land_tint_color']]
        water_tint_color  = [int(c * 255) for c in kwargs['water_tint_color']]

        ax      = Axes(self.ds)
        xsize   = kwargs['geometry'][0]
        ysize   = kwargs['geometry'][1]
        xinches = float(ax.xhigh - ax.xlow)
        yinches = float(ax.yhigh - ax.ylow)
        xpixels = int(round(xinches/ax.xsize * xsize))
        ypixels = int(round(yinches/ax.ysize * ysize))
        xpos    = int(round(ax.xlow / ax.xsize * xsize))
        ypos    = int(round((ax.ysize - ax.yhigh) / ax.ysize * ysize))

#       Create the map image

        buf = io.BytesIO()
        plt.gca().set_axis_off()
        plt.savefig(buf, format='png', bbox_inches='tight', pad_inches=0,
                        dpi=300, facecolor=face_color)

        buf.seek(0)
        fg = Image.open(buf)

        land  = self.immask(fg,0)
        water = self.immask(fg,255)

        land  = self.imfill(land, 0, land_color)
        land  = self.imtint(land, land_tint_color)
        land  = self.imbright(land, land_brightness)

        water = self.imfill(water, 255, water_color)
        water = self.imtint(water, water_tint_color)
        water = self.imbright(water, water_brightness)

        fg.paste(land, (0,0), land)
        fg.paste(water, (0,0), water)

#       Re-size and position the map image on a background image
#       matching the size of the final GrADS image.

        fg = fg.resize((xpixels, ypixels), Image.ANTIALIAS)

        if isGrayscale: 
            bg = Image.new("RGB", [xsize, ysize], face_color).convert('LA')
        else:
            bg = Image.new("RGB", [xsize, ysize], face_color)

        bg.paste(fg, (xpos, ypos), fg)
        bg.save(fname, "PNG")
        bg.close()

    def immask(self, src, value):

        data    = src.getdata()
        mask    = self.mask.getdata()
        img     = src.copy()
        newData = []

        for i in range(0,len(data)):
            idata = data[i]
            imask = mask[i]

            if imask[0] == value:
                newData.append(idata)
            else:
                newData.append((0,0,0,0))

        img.putdata(newData)

        return img

    def imfill(self, src, value, color):

        if self.is_clear(color): return src

        data    = src.getdata()
        mask    = self.mask.getdata()
        img     = src.copy()
        newData = []
        color   = tuple(color)

        for i in range(0,len(data)):
            idata = data[i]
            imask = mask[i]

            if imask[0] == value:
                newData.append(color)
            else:
                newData.append(idata)

        img.putdata(newData)

        return img

    def imtint(self, src, tint_color):

        if self.is_clear(tint_color): return src

        tint_color  = '#%02x%02x%02x'%tuple(tint_color)
        return self.image_tint(src, tint_color)

    def imbright(self, src, brightness):

        if brightness == 1.0: return src

        enhancer = ImageEnhance.Brightness(src)
        return enhancer.enhance(brightness)

    def is_solid(self, color): 
        if len(color) < 4: return True
        if color[3] > 0.0: return True
        return False

    def is_clear(self, color):
        return not self.is_solid(color)

    def image_tint(self, src, tint='#ffffff'):

        if Image.isStringType(src):  # file path?
            src = Image.open(src)
        if src.mode not in ['RGB', 'RGBA']:
            raise TypeError('Unsupported source image mode: {}'.format(src.mode))
        src.load()

        tr, tg, tb = getrgb(tint)
        tl = getcolor(tint, "L")  # tint color's overall luminosity
        if not tl: tl = 1  # avoid division by zero
        tl = float(tl)  # compute luminosity preserving tint factors
        sr, sg, sb = map(lambda tv: tv/tl, (tr, tg, tb))  # per component
                                                      # adjustments
        # create look-up tables to map luminosity to adjusted tint
        # (using floating-point math only to compute table)
        luts = (tuple(map(lambda lr: int(lr*sr + 0.5), range(256))) +
                tuple(map(lambda lg: int(lg*sg + 0.5), range(256))) +
                tuple(map(lambda lb: int(lb*sb + 0.5), range(256))))
        l = grayscale(src)  # 8-bit luminosity version of whole image
        if Image.getmodebands(src.mode) < 4:
            merge_args = (src.mode, (l, l, l))  # for RGB verion of grayscale
        else:  # include copy of src image's alpha layer
            a = Image.new("L", src.size)
            a.putdata(src.getdata(3))
            merge_args = (src.mode, (l, l, l, a))  # for RGBA verion of grayscale
            luts += tuple(range(256))  # for 1:1 mapping of copied alpha values

        return Image.merge(*merge_args).point(luts)

    def add_map_layer(self, map, bmap):

        service     = bmap['service']

        if service == 'shaderelief':
            map.shadedrelief()
        elif service == 'bluemarble':
            map.bluemarble()
        elif service == 'etopo':
            map.etopo()
        elif service == 'drawlsmask':
            print 'Using cached mask'
        else:
            map.arcgisimage(service=service, xpixels = 1000)

    def create_mask(self, map, **kwargs):

#       Create the land/sea mask

        map.drawlsmask(land_color='black',ocean_color='white',
                            resolution='f', grid=1.25, lakes=True)

        buf = io.BytesIO()
        plt.savefig(buf, format='png', bbox_inches='tight', pad_inches=0,
                        dpi=300, facecolor='red')

        buf.seek(0)
        self.mask = Image.open(buf).convert('LA')

        plt.clf()

    def get_map_handle(self, **kwargs):

        mproj      = kwargs['mproj']
        mpvals     = kwargs['mpvals']
        xsize      = kwargs['geometry'][0]
        ysize      = kwargs['geometry'][1]
        lat1, lat2 = kwargs['lat']
        lon1, lon2 = kwargs['lon']
        if len(mpvals) == 4: lon1, lon2, lat1, lat2 = mpvals

        lon_center = lon1 + 180

        plt.clf()
        plt.figure(figsize=(xsize*2/300.0, ysize*2/300.0), dpi=100)
        plt.gca().set_axis_off()

        if mproj == 'nps' and lat2 == 90:

            map = Basemap(resolution=None, projection='npstere',
                          boundinglat=lat1,lon_0=lon_center, round=True)

        elif mproj == 'nps':

            map = self.polar_stere(lon1, lon2, lat1, lat2)

        elif mproj == 'sps' and lat1 == -90:

            map = Basemap(resolution=None, projection='spstere',
                          boundinglat=lat2,lon_0=180, round=True)

        elif mproj == 'sps':

            map = self.polar_stere(lon1, lon2, lat1, lat2)

        elif mproj == 'orthogr':

            map = Basemap(resolution=None, projection='ortho',
                    lon_0=(lon1+lon2)/2.0, lat_0=0.0, round=True)

        else:

            if lon1 > 180.:
                lon1 -= 360.
                lon2 -= 360.

            map = Basemap(resolution=None, projection='cyl',
                          llcrnrlat=lat1,urcrnrlat=lat2,
                          llcrnrlon=lon1,urcrnrlon=lon2)

        return map

    def polar_stere(self, lon_w, lon_e, lat_s, lat_n, **kwargs):
        '''Returns a Basemap object (NPS/SPS) focused in a region.

        lon_w, lon_e, lat_s, lat_n -- Graphic limits in geographical coordinates
                                      W and S directions are negative.
        **kwargs -- Aditional arguments for Basemap object.


        '''

#       Code acquired from:
#       http://code.activestate.com/recipes/

        lon_0 = lon_w + (lon_e - lon_w) / 2.
        ref = lat_s if abs(lat_s) > abs(lat_n) else lat_n
        lat_0 = math.copysign(90., ref)
        proj = 'npstere' if lat_0 > 0 else 'spstere'
        prj = Basemap(projection=proj, lon_0=lon_0, lat_0=lat_0,
                              boundinglat=0, resolution='c')
        #prj = pyproj.Proj(proj='stere', lon_0=lon_0, lat_0=lat_0)
        lons = [lon_w, lon_e, lon_w, lon_e, lon_0, lon_0]
        lats = [lat_s, lat_s, lat_n, lat_n, lat_s, lat_n]
        x, y = prj(lons, lats)
        ll_lon, ll_lat = prj(min(x), min(y), inverse=True)
        ur_lon, ur_lat = prj(max(x), max(y), inverse=True)
        return Basemap(projection='stere', lat_0=lat_0, lon_0=lon_0,
                               llcrnrlon=ll_lon, llcrnrlat=ll_lat,
                               urcrnrlon=ur_lon, urcrnrlat=ur_lat, **kwargs)

    def get_map_key(self, bmaps):

        keystr = ''
        basemap_keys = ['lat', 'lon', 'mpvals', 'mproj',
                        'geometry', 'lights_off', 'grayscale',
                        'land_color', 'water_color', 'ocean_color',
                        'lake_color', 'tint_color', 'land_tint_color',
                        'water_tint_color', 'brightness', 'service',
                        'land_brightness', 'water_brightness']

        for bmap in bmaps:
            d = {k:bmap[k] for k in basemap_keys}
            keystr += json.dumps(d)

        keystr = keystr.replace("u'", "'")
        return str(uuid.uuid3(uuid.NAMESPACE_DNS,keystr))

    def make_labels(self):

        self.ds('query gxinfo')

        xlo = float(self.ds.rword(3,4))
        xhi = float(self.ds.rword(3,6))
        ylo = float(self.ds.rword(4,4))
        yhi = float(self.ds.rword(4,6))

        bbox = (xlo, ylo, xhi, yhi)

        self.draw_labels(bbox,  1, self.lang.header_labels)
        self.draw_labels(bbox, -1, self.lang.trailer_labels)
        self.draw_labels(bbox,  1, self.lang.special_labels)

    def draw_labels(self, bbox, order, labels):

        if not labels: return

        factor = 1
        ypos   = bbox[3]
        
        if order == -1:
            factor = -1
            ypos   = bbox[1]

        xmid = (bbox[2] - bbox[0]) / 2.0 + bbox[0]

        for type in labels:

            cmds = self.labels.get(type, [])
            if not cmds: continue
          
            attr    = [ c.split() for c in cmds if self.lang.is_attribute(c)]
            attr    = { c[0][1:]:c[1:] for c in attr }
            cmds    = [ c.split(' ') for c in cmds if self.is_valid(c) ]
            params  = { c[1]:c[2:] for c in cmds if c[0] == 'set' }
            drawcmd = [ c for c in cmds if c[0] == 'draw' ]
            cmds    = [ ' '.join(c) for c in cmds if c[0] != 'draw' ]
            cmds    = '\n'.join(cmds)

            self.ds(cmds)

            drawcmd = drawcmd[-1]
            size    = float(params['strsiz'][-1])
            pos     = params['string'][1]
            margin  = float(attr['MARGIN'][0])
            string  = ' '.join(drawcmd[2:]).strip('\\')

            xpos = xmid
            if pos == 'l': xpos = bbox[0]
            if pos == 'r': xpos = bbox[2]
            ypos = ypos + (margin + 0.5 * size) * factor

            self.ds('draw string ' + str(xpos) + ' ' + str(ypos) + ' ' + string)

            ypos = ypos + (0.5 * size) * factor

    def draw_subtitle(self, obj):

        cmds = [ cmd for cmd in obj.cmds[0:-1] if self.is_valid(cmd) ]
        cmds = '\n'.join(cmds)
        self.ds(cmds)

        cmd   = re.sub("\s+"," ",obj.cmds[-1].strip()).split()
        title = ' '.join(cmd[2:])
        title = title.split('|')

        self.ds('query gxinfo')

        xlo = float(self.ds.rword(3,4))
        xhi = float(self.ds.rword(3,6))
        ylo = float(self.ds.rword(4,4))
        yhi = float(self.ds.rword(4,6))

        xmid = (xhi - xlo) / 2.0 + xlo
        ymid = yhi + 0.6

        self.ds('draw string ' + str(xmid) + ' ' + str(ymid) + ' ' + title[0])

        if len(title) <= 1: return

        ymid = ymid - 0.20
        self.ds('draw string ' + str(xmid) + ' ' + str(ymid) + ' ' + title[1])

    def draw_tmstring(self, obj):

        cmds = [ cmd for cmd in obj.cmds[0:-1] if self.is_valid(cmd) ]
        cmds = '\n'.join(cmds)
        self.ds(cmds)

        cmd   = obj.cmds[-1].strip()
        title = cmd[13:]

        self.ds('query gxinfo')

        xlo = float(self.ds.rword(3,4))
        xhi = float(self.ds.rword(3,6))
        ylo = float(self.ds.rword(4,4))
        yhi = float(self.ds.rword(4,6))

        xmid = (xhi - xlo) / 2.0 + xlo
        ymid = yhi + 0.17

        self.ds('draw string ' + str(xmid) + ' ' + str(ymid) + ' ' + title)

    def w2xy(self, wx, wy):

        self.ds('query w2xy %s %s'%(wx, wy),Quiet=True)
        x = float(self.ds.rword(1,3))
        y = float(self.ds.rword(1,6))

        return (x, y)

    def is_valid(self, cmd):

        if self.lang.is_define(cmd) or \
               self.lang.is_macro(cmd) or \
                   self.lang.is_attribute(cmd) or \
                             self.lang.is_end(cmd) or \
                                 self.lang.is_addon(cmd):

            return False

        return True

    def get_state(self, state, key, default=novalue):

        if key not in state:
            if default is novalue: return ''
            return default

        return state[key]

class Axes(object):

    def __init__(self, dataservice, **kwargs):

        self.ds      = dataservice

        self.ds('query gxinfo')

        xsize = float(self.ds.rword(2,4))
        ysize = float(self.ds.rword(2,6))
        ylow  = float(self.ds.rword(4,4))
        yhi   = float(self.ds.rword(4,6))
        xlow  = float(self.ds.rword(3,4))
        xhi   = float(self.ds.rword(3,6))

        self.LL = (0.0, 0.0);   self.Ll = (xlow, 0.0);   self.Lr = (xhi, 0.0)
        self.LR = (xsize, 0.0); self.lL = (0.0, ylow);   self.ll = (xlow, ylow)
        self.lr = (xhi, ylow);  self.lR = (xsize, ylow); self.uL = (0.0, yhi)
        self.ul = (xlow, yhi);  self.ur = (xhi, yhi);    self.uR = (xsize, yhi)
        self.UL = (0.0, ysize); self.Ul = (xlow, ysize); self.Ur = (xhi, ysize)
        self.UR = (xsize, ysize)

        self.xsize = xsize
        self.ysize = ysize
        self.ylow  = ylow
        self.yhigh = yhi
        self.xlow  = xlow
        self.xhigh = xhi

        self.clipit  = kwargs.get('clip', True)

    def get(self, key, default=None): return self.__dict__.get(key, default)

    def clip(self, line, p=None):

        if not self.clipit: return line

        p1   = line[0]
        p2   = line[1]
        if p is None: return self.clip(line, p1) + self.clip(line, p2)

        if not self.is_clipped(p): return [p]

      # Clip the point to the nearest
      # axis boundary.

        xb = max(p[0], self.xlow)
        xb = min(xb, self.xhigh)
        yb = max(p[1], self.ylow)
        yb = min(yb, self.yhigh)

      # Simple check if slope is zero
      # or undefined.

        if p1[0] == p2[0] or p1[1] == p2[1]:
            if self.is_extrapolated(line, (xb,yb)): return []
            return [(xb,yb)]

      # Determine new coordinates based
      # on where the line intersects the
      # axes.

        slope = (p2[1] - p1[1]) / (p2[0] - p1[0])
        intercept = p1[1] - slope*p1[0]

        p = (xb, slope*xb + intercept)
        if not self.is_clipped(p):
            if self.is_extrapolated(line, p): return []
            return [p]

        p = ((yb - intercept) / slope, yb)
        if not self.is_clipped(p):
            if self.is_extrapolated(line, p): return []
            return [p]

        return []

    def is_clipped(self, p):

        if not self.clipit: return False

        if p[0] < self.xlow: return True
        if p[0] > self.xhigh: return True
        if p[1] < self.ylow: return True
        if p[1] > self.yhigh: return True

        return False

    def is_extrapolated(self, line, p):

        p1   = line[0]
        p2   = line[1]

        extrapolated = True
        if p[0] >= p1[0] and p[0] <= p2[0]: extrapolated = False
        if p[0] <= p1[0] and p[0] >= p2[0]: extrapolated = False

        if extrapolated: return True
        
        if p[1] >= p1[1] and p[1] <= p2[1]: return False
        if p[1] <= p1[1] and p[1] >= p2[1]: return False

        return True

class Colorbar(object):

    def __init__(self, dataservice):

        self.cbars = []

        self.ds = dataservice
        self.ax = Axes(self.ds)

        self.default = DictContainer()

        self.default.str_offset      = (0.05, 0.05)
        self.default.tick_label_size = (0.12, 0.13)
#       self.default.tick_label_size = (0.20, 0.21)
        self.default.bar_label_size  = (0.12, 0.12)
#       self.default.bar_label_size  = (0.20, 0.20)
        self.default.bar_offset      = 0.075
        self.default.bar_thick       = 0.15

    def add(self, colors, levels, **kwargs):

        cbar = DictContainer()
        cbar.update(self.default)

        cbar.colors  = list(colors)
        cbar.levels  = list(levels)
        cbar.options = dict(kwargs)

        self.cbars.append(cbar)

class HorizontalColorbar(Colorbar):

    def draw(self, **kwargs):

        ax = self.ax
        self.set_bar_dims(**kwargs)
        self.set_label_sizes(**kwargs)

        for index, cbar in enumerate(self.cbars):

            colors  = cbar.colors
            ncolors = len(colors)
            x, y    = cbar.beg

            for n in range(ncolors):

                xwid = cbar.stride[0]
                if n == 0 or n == ncolors-1: xwid = xwid / 2.0

                xl = x; xr = x + xwid; yb = y; yt = y + cbar.size[1]

                self.ds('set line %s'%(colors[n]))
                self.ds('draw recf %f %f %f %f'%(xl,yb,xr,yt))

                x += xwid

            self.ds('set line 1 1 5')
            self.ds('draw rec %f %f %f %f'%(cbar.beg+cbar.end))

            self.draw_tick_labels(cbar)

            label = cbar.options.get('cblabel',None)

            xp = cbar.beg[0]
            yp = cbar.end[1] + cbar.str_offset[1]
            if label:
                self.ds('set font 13')
                self.ds('set strsiz %f %f'%cbar.bar_label_size)
                self.ds('set string 1 bl 5 0')
                self.ds('draw string %f %f %s'%(xp,yp,label))
                self.ds('set string 1 bl 5 0')

    def set_bar_dims(self, **kwargs):

        tight     = kwargs.get('tight', False)
        cbar_only = kwargs.get('cbar_only', False)

        for index, cbar in enumerate(self.cbars):

            ax      = self.ax
            nbars   = len(self.cbars)
            ncolors = len(cbar.colors)

            length = ax.xhigh - ax.xlow
            if cbar_only: length = 10
            offset = (1.0 - 0.9) * (length / nbars) / 2.0
            xbeg   = (length / nbars) * index + ax.xlow + offset
            ybeg   = ax.ylow - 0.9
            if tight: ybeg = ax.ylow - 0.5
            if cbar_only: ybeg = (ax.yhigh - ax.ylow) / 2.0
            stride = (length / nbars - offset * 2) / (ncolors - 1)
            xend   = xbeg + stride * (ncolors - 1)
            yend   = ybeg + self.default.bar_thick

            cbar.size   = (xend-xbeg, self.default.bar_thick)
            cbar.beg    = (xbeg, ybeg)
            cbar.end    = (xend, yend)
            cbar.offset = (offset, 0)
            cbar.stride = (stride, 0)

    def set_label_sizes(self, **kwargs):

        scale         = float(kwargs.get('scale', 1.0))
        min_char_size = [sz*scale for sz in self.default.tick_label_size]

        for cbar in self.cbars:

            mxlen   = 0
            nlevels = len(cbar.levels)
            skip    = max(int(cbar.options.get('skip',1)),1)

            for n in range(0, nlevels, skip):

                label = '%0.4f'%(float(cbar.levels[n]))
                label = label.strip('0').rstrip('.')
                if label == '': label = '0'
                mxlen  = max(mxlen, len(label))

            char_size = cbar.stride[0]*skip / float(mxlen)
            char_size = min(char_size, self.default.tick_label_size[0])
            char_size = char_size * scale
            cbar.tick_label_size = (char_size, 1.1*char_size)

            min_char_size = min(min_char_size, char_size)

        for cbar in self.cbars:
            char_size = float(cbar.options.get('tklabsiz', min_char_size))
            cbar.tick_label_size = (char_size, 1.1*char_size)
            cbar.bar_label_size  = [sz*scale for sz in cbar.bar_label_size]
            cbar.bar_label_size  = tuple(cbar.bar_label_size)

    def draw_tick_labels(self, cbar):

        x, y    = cbar.beg
        levels  = cbar.levels
        nlevels = len(levels)
        skip    = max(int(cbar.options.get('skip',1)),1)

        for n in range(nlevels):

            xwid = cbar.stride[0]
            if n == 0: xwid = xwid / 2.0

            xl = x; xr = x + xwid; yb = y; yt = y + cbar.size[1]
            yp  = yb - cbar.str_offset[1]
            ytk = (yt-yb) / 3.0 + yb

            if n%skip == 0:

                label = '%0.4f'%(float(levels[n]))
                label = label.strip('0').rstrip('.')
                if label == '': label = '0'
                if label == '-0': label = '0'

                self.ds('set font 13')
                self.ds('set strsiz %f %f'%cbar.tick_label_size)
                self.ds('set string 1 tc 5 0')
                self.ds('draw string %f %f %s'%(xr,yp,label))

                self.ds('set line 1 1 5')
#               self.ds('draw line %f %f %f %f'%(xl,yb,xl,ytk))
                self.ds('draw line %f %f %f %f'%(xr,yb,xr,ytk))

            x += xwid

class VerticalColorbar(Colorbar):

    def draw(self, **kwargs):

        ax = self.ax
        self.set_bar_dims(**kwargs)
        self.set_label_sizes(**kwargs)

        for index, cbar in enumerate(self.cbars):

            colors  = cbar.colors
            ncolors = len(colors)
            x, y    = cbar.beg

            for n in range(ncolors):

                stride = cbar.stride[1]
                if n == 0 or n == ncolors-1: stride = stride / 2.0

                yb = y; yt = y + stride; xl = x; xr = x + cbar.size[0]

                self.ds('set line %s'%(colors[n]))
                self.ds('draw recf %f %f %f %f'%(xl,yb,xr,yt))

                y += stride

            self.ds('set line 1 1 5')
            self.ds('draw rec %f %f %f %f'%(cbar.beg+cbar.end))

            self.draw_tick_labels(cbar)

            label = cbar.options.get('cblabel',None)

            xp = cbar.beg[0] - cbar.str_offset[0]
            yp = cbar.beg[1] + (cbar.end[1] - cbar.beg[1]) / 2.0
            if label:
                self.ds('set font 13')
                self.ds('set strsiz %f %f'%cbar.bar_label_size)
                self.ds('set string 1 bc 5 90')
                self.ds('draw string %f %f %s'%(xp,yp,label))
                self.ds('set string 1 bc 5 0')

    def set_bar_dims(self, **kwargs):

        for index, cbar in enumerate(self.cbars):

            ax        = self.ax
            nbars     = len(self.cbars)
            ncolors   = len(cbar.colors)
            bar_label = cbar.options.get('cblabel', None)

            length = ax.yhigh - ax.ylow
            offset = (1.0 - 0.9) * (length / nbars) / 2.0

            xbeg   = ax.xhigh + 0.1
            if bar_label: xbeg += cbar.bar_label_size[0] + cbar.str_offset[0]

            ybeg   = (length / nbars) * index + ax.ylow + offset
            stride = (length / nbars - offset * 2) / (ncolors - 1)
            xend   = xbeg + self.default.bar_thick
            yend   = ybeg + stride * (ncolors - 1)

            cbar.size   = (self.default.bar_thick, yend-ybeg)
            cbar.beg    = (xbeg, ybeg)
            cbar.end    = (xend, yend)
            cbar.offset = (0, offset)
            cbar.stride = (0, stride)

    def set_label_sizes(self, **kwargs):

        ax = self.ax

        scale         = float(kwargs.get('scale', 1.0))
        min_char_size = [sz*scale for sz in self.default.tick_label_size]

        for cbar in self.cbars:

            mxlen   = 0
            nlevels = len(cbar.levels)
            margin  = ax.xsize - cbar.end[0]
            length  = cbar.size[1]
            skip    = max(int(cbar.options.get('skip',1)),1)

            for n in range(0, nlevels, skip):

                label = '%0.4f'%(float(cbar.levels[n]))
                label = label.strip('0').rstrip('.')
                if label == '': label = '0'
                mxlen = max(mxlen, len(label))

            nlevels   = round(float(nlevels) / skip)
            char_size = margin / float(mxlen)
            char_size = min(char_size, self.default.tick_label_size[0])
            char_size = min(char_size, length / float(nlevels))
            char_size = char_size * scale
            cbar.tick_label_size = (char_size, char_size)

            min_char_size = min(min_char_size, char_size)

        for cbar in self.cbars:
            char_size = float(cbar.options.get('tklabsiz', min_char_size))
            cbar.tick_label_size = (char_size, char_size)
            cbar.bar_label_size  = [sz*scale for sz in cbar.bar_label_size]
            cbar.bar_label_size  = tuple(cbar.bar_label_size)

    def draw_tick_labels(self, cbar):

        x, y    = cbar.beg
        levels  = cbar.levels
        nlevels = len(levels)
        skip    = max(int(cbar.options.get('skip',1)),1)

        for n in range(nlevels):

            stride = cbar.stride[1]
            if n == 0: stride = stride / 2.0

            yb = y; yt = y + stride; xl = x; xr = x + cbar.size[0]
            xp  = xr + cbar.str_offset[0]
            xtk = xr - cbar.size[0] / 3.0

            if n%skip == 0:

                label = '%0.4f'%(float(levels[n]))
                label = label.strip('0').rstrip('.')
                if label == '': label = '0'

                self.ds('set font 13')
                self.ds('set strsiz %f %f'%cbar.tick_label_size)
                self.ds('set string 1 l 5 0')
                self.ds('draw string %f %f %s'%(xp,yt,label))

                self.ds('set line 1 1 5')
                self.ds('draw line %f %f %f %f'%(xtk,yt,xr,yt))


            y += stride

class DictContainer(object):

    def __init__(self, **kwargs):

        self.__dict__.update(kwargs)

    def update(self, defs):

        if isinstance(defs, DictContainer):
            self.__dict__.update(defs.__dict__)
            return

        if isinstance(defs, dict):
            self.__dict__.update(defs)
            return
