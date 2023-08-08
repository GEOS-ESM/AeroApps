#--------------------------------------------------------------------------
#
#  REVISION HISTORY:
#
#  18Dec2007  da Silva  First crack.
#
#--------------------------------------------------------------------------
#
#    Copyright (C) 2007-2008 by Arlindo da Silva <dasilva@opengrads.org>
#    All Rights Reserved.
#
#    This program is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation# using version 2 of the License.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY# without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program# if not, please consult  
#              
#              http://www.gnu.org/licenses/licenses.html
#
#    or write to the Free Software Foundation, Inc., 59 Temple Place,
#    Suite 330, Boston, MA 02111-1307 USA
#
#------------------------------------------------------------------------

"""
This module implements the GaLab class, an extention of the GrADS
client class which adds Matplotlib/Basemap specific methods for
contours, images and other graphical functionality. This class
provides high level methods operating directly on GrADS expressions
(or *fields*) while retaining all the configurability that Matplotlib
has to offer.
"""

__version__ = '1.2.2'

from gacore    import dt2gat, gat2dt
from ganum     import *
from numtypes  import *
from tempfile  import mktemp
from datetime  import datetime, timedelta

import matplotlib

from matplotlib.image import imread

from numpy             import meshgrid, arange, shape, size, reshape, tile
from matplotlib.pyplot import show, draw, colorbar, axes, axis, gcf, cm, \
                              imshow, get_cmap, figimage, quiverkey, \
                              isinteractive

try:
    from mpl_toolkits.basemap import Basemap, interp
except:
    from matplotlib.toolkits.basemap import Basemap, interp

try:
    from PIL              import Image
    from matplotlib.image import pil_to_array
    HAS_PIL = True
except ImportError:
    HAS_PIL = False

def _pil_to_array(im):
    """
    Wrapper to handle mod introduced in matplotlib 1.2 that turned array upside down.
    """
    rgb = pil_to_array(im)
    M, m, r = matplotlib.__version__.split('.')
    ver = int(M)*100 + int(m)
    if ver>101:  # matplotlib 1.2 and later
        rank = len(rgb.shape)
        if rank == 2:
            rgb = rgb[-1::-1,:]
        else:
            rgb = rgb[-1::-1,:,:]
    return rgb
            
class GaLab(GaNum):
    """
    This class extends the GrADS client class GaNum by providing an interface
    to PyLab and Matplotlib. The main design philosophy is to use the
    GrADS dimension environment to derive sensible defaults for the
    high level graphical commands, while retaining most of Matplotlib
    configurability.
    """

    def __init__ (self, **kwopts):
        """ 
        Contructor. All keyword options (**kwopts) are passed to the
        costructor of the core GrADS client class.
        """
  
        GrADS.__init__(self, **kwopts)

        self.map = None         # basemap
        self.map_stamp = None

        self.blue = None         # blue marble image 
        self.blue_stamp = None
        self.blue_size = (-1,-1)

        self.datadir = os.sep.join([os.path.dirname(__file__), 'data'])
        self.blue_file = os.sep.join([self.datadir, 'blue_marble.jpg'])

        self.transf = True # usually do transforms
        self.MAPRES = 'l'  # Map default resolution
        self.AREA_THRESH = 10000.
#   ..................................................................

    def basemap ( self, proj='latlon',dh=None,opts=None,
                  resolution='', area_thresh='', **kwopts ):
        """
        Defines basemap based on GrADS dimension environment.
        Both *resolution* and *area_thresh* have sticking defaults,
        starting with 'l' and 10000.
        """

        # Sticking defaults
        # -----------------
        if resolution == '':
            resolution = self.MAPRES
        else:
            self.MAPRES = resolution # it sticks
        if area_thresh == '':
            area_thresh = self.AREA_THRESH
        else:
            self.AREA_THRESH = area_thresh 
 
        if dh==None:
            dh = self.query("dims", Quiet=True)

#       Coping with roundoff...
#       -----------------------
        dh.lat = list(dh.lat)
        if dh.lat[0] < -90.: dh.lat[0] = -90.
        if dh.lat[1] >  90.: dh.lat[1] =  90.

        self.axis = [0.1,0.1,0.75,0.75]

        if proj=='latlon' or proj=='cyl':
            self.proj = 'cyl'
            self.labels = [1,0,0,1]
            if opts == None:
                self.transf = False
                self.dh = dh
                ll_lon = dh.lon[0]
                ll_lat = dh.lat[0]
                ur_lon = dh.lon[1]
                ur_lat = dh.lat[1]
            elif len(opts)==4:
                self.transf = True
                ll_lon,ll_lat,ur_lon,ur_lat = opts
            else:
                raise GrADSError, \
                'Projection "geos" expects opts = ' + \
                ' (ll_lon,ll_lat,ur_lon,ur_lat,lon_0,satellite_height)'
            self.map = Basemap(llcrnrlon=ll_lon, llcrnrlat=ll_lat, \
                                urcrnrlon=ur_lon, urcrnrlat=ur_lat, \
                                resolution=resolution, area_thresh=area_thresh,\
                                projection='cyl', **kwopts)
        elif proj=='geos' or proj=='geos_':
            self.proj = proj
            self.transf = True
            self.labels = [0,0,0,0]
            if opts == None:
                ll_lon = None
                ll_lat = None
                ur_lon = None
                ur_lat = None
                lon_0 = -75. # subpoint longitude (GOES west)
                
            elif (type(opts) is IntType) or (type(opts) is FloatType) :
                ll_lon = None
                ll_lat = None
                ur_lon = None
                ur_lat = None
                lon_0 = opts  # subpoint longitude

            elif len(opts)==5:
                ll_lon,ll_lat,ur_lon,ur_lat,lon_0 = opts

            else:
                raise GrADSError, \
                'Projection "geos" expects opts=subpoint_lon  or' + \
                ' opts=(ll_lon,ll_lat,ur_lon,ur_lat,subpoint_lon)'

            if proj=='geos':
                self.map = Basemap(projection='geos', \
                                llcrnrlon=ll_lon, llcrnrlat=ll_lat, \
                                urcrnrlon=ur_lon, urcrnrlat=ur_lat, 
                                lon_0=lon_0, \
                                resolution=resolution, area_thresh=area_thresh, \
                                rsphere=(6378137.00,6356752.3142),
                                satellite_height = 35785831.0,
                                **kwopts )
            else:
                self.map = Basemap(projection='geos', \
                                resolution=resolution, area_thresh=area_thresh, \
                                rsphere=(6378137.00,6356752.3142),
                                satellite_height = 35785831.0,
                                **kwopts )
                

                
        elif proj=='orthogr' or proj=='ortho' or proj=='npo' or proj=='spo': 
            self.proj = 'ortho'
            self.transf = True
            self.labels = [0,0,0,0] 
            if opts == None:
                if proj=='npo':
                    lon_0 = -100.
                    lat_0 = +90.
                elif proj=='spo':
                    lon_0 = -50.
                    lat_0 = -90.
                elif dh.lat[0] >= 0:
                    lon_0 = -100.
                    lat_0 = +90.
                else:
                    lon_0 = -90.
                    lat_0 = 45.
            elif len(opts)==2:
                lon_0 = opts[0]
                lat_0 = opts[1]
            else:
                raise GrADSError, \
                'Projection "ortho" expects opts = (lon_0,lat_0)'

            self.map = Basemap(projection='ortho', \
                                lon_0=lon_0, lat_0=lat_0, \
                                resolution=resolution, area_thresh=area_thresh,\
                                **kwopts )

        elif proj=='stereo' or proj=='stere' or proj=='nps' or proj=='sps':
            self.proj = 'stereo'
            self.transf = True
            self.labels = [0,0,0,0]
            self.axis = [0.05,0.1,0.8,0.8]
            if opts == None:
                if proj=='nps':
                    lon_0 = -90.
                    lat_0 = +20.
                elif proj=='sps':
                    lon_0 = -50.
                    lat_0 = -20.
                elif dh.lat[0] >= 0:
                    lon_0 = -90.
                    lat_0 = +20.
                else:
                    lon_0 = -50.
                    lat_0 = -20.
            elif len(opts)==2:
                lon_0 = opts[0]
                lat_0 = opts[1]
            else:
                raise GrADSError, \
                'Projection "%s" expects opts = (lon_0,boundinglat)'%proj

            if lat_0>0:
                proj1 = 'npstere'
            else:
                proj1 = 'spstere'
                
            self.map = Basemap(projection=proj1, \
                                lon_0=lon_0, boundinglat=lat_0, \
                                resolution=resolution, area_thresh=area_thresh,\
                                **kwopts )

        else:

            self.proj = proj
            self.transf = True
            self.labels = [0,0,0,0]
            self.map = Basemap(projection=proj, \
                                resolution=resolution, area_thresh=area_thresh, \
                                **kwopts )

        self.map_stamp = time()  # basemap timestamp
        
#   ..................................................................

    def blue_marble(self, mode=None, Filename=None, dh=None, 
                          Show=False, Nx=None, Ny=None, Cache=True):
        """ 
            Loads Blue Marble image at the current map projection and
            use it as the default background image. For turning the
            Blue Marble background on for subsequent plotting methods
            you enter

               ga.blue_marble('on')

            For retriving an image on the current map projection:

               bm = blue_marble()

        """

        if not HAS_PIL:
            raise GrADSError, 'PIL is required to for method "blue_marble"'

#       User specified file name
#       ------------------------
        new_file = False
        if (Filename is not None) and (Filename is not self.blue_file):
            self.blue_file = Filename
            new_file = True

#       Set the mode on/off
#       -------------------
        if mode is not None:
	    if mode=='off' or mode=='OFF' or mode=='Off':
                self.blue = None
 	    elif mode=='on' or mode=='ON' or mode=='On':
	        self.blue = ma.zeros((1,1,4)) # Will be filled in at first use
	    else:
                raise GrADSError, 'Invalid Blue Marble mode'
	    self.blue_stamp = None # will force update on first use
	    self.blue_size = (-1,-1)            
            return None

#       Make sure we have a basemap
#       ---------------------------
        if self.map == None or self.transf == False:
            self.basemap(dh=dh)
        m = self.map

#       Default size based on map projection
#       ------------------------------------
        if Ny==None:
            Ny = 600
        if Nx==None:
            aspect = (m.xmax-m.xmin)/(m.ymax-m.ymin)
            Nx = max(450,int(aspect * Ny))

#       If the right stuff is in cache, return it
#       -----------------------------------------
        if not new_file:
            if self.blue_stamp==self.map_stamp and self.blue_size==(Nx,Ny):
                if Show:
                    self.map.imshow(self.blue)
                    return None
                else:
                    return self.blue

#       Read image and make an RGBA array out of it
#       -------------------------------------------
        pim = Image.open(self.blue_file)
        if pim.mode is 'RGB':
            try:
                pim.putalpha(255)
            except:
                pass
        rgba = _pil_to_array(pim).astype(float32)/255. # normalized floats.


#       Define lat/lon grid that image spans (projection='cyl').
#       --------------------------------------------------------
        nlons = rgba.shape[1]; nlats = rgba.shape[0]
        delta = 360./float(nlons)
        dh = delta / 2.
        lons = arange(-180.+dh,180.+dh,delta)
        lats = arange(-90.+dh,90+dh,delta)

#       Shift image longitudinally to match projection (cyl-like only)
#       -------------------------------------------------------------
        if self.proj in ( 'cyl', 'merc', 'mill'):
            if self.map.xmin != -180.:
                d = self.map.xmin + 180.
                lons = lons + d 
                i = int(d/delta)
                n = nlons - i 
                if i<0 or i>=nlons:
                    raise GrADSError, "Internal error shifting image, i=%d"%i
                tmp = rgba.copy()
                rgba[:,0:n,:] = tmp[:,i:nlons,:]
                rgba[:,n:nlons,:] = tmp[:,0:i,:]
                
#       Transform the data to current map transformation
#       ------------------------------------------------
        bm = ma.zeros((Ny,Nx,4))
        for k in range(4):
            bm[:,:,k] = m.transform_scalar(rgba[:,:,k],lons,lats,Nx,Ny,
                                           masked=True)
        bm = bm.filled(fill_value=0.) # no longer masked

#       Cache it
#       --------
        if Cache:
            self.blue = bm
            self.blue_stamp = self.map_stamp 
            self.blue_size = (Nx,Ny)

#       All done
#       --------
        if Show:
            self.map.imshow(bm)
            return None
        else:
            return bm

#   ..................................................................

    def implain(self, expr, dh=None, cmap=None, interpolation='bicubic', \
               vmin=None, vmax=None, **kwopts):
        """Evaluates GrADS expression and display it as an image
        on the current basemap."""

#       Check dim environment
#       ---------------------
        if dh==None:
            dh = self.query("dims", Quiet=True)
        if dh.nz>1 or dh.nt>1:
            raise GrADSError, 'Not a horizontal slice; ' + \
                  'expected (nz,nt)=(1,1) but got (%d,%d)'%(dh.nz,dh.nt)

#       Evaluate GrADS expression
#       -------------------------
        X = self.exp(expr)
        g = X.grid

#       When vmin/vmax are specified, make the over/under colors
#       transparent
#       -------------------------------------------------------
        if vmin!=None and vmax!=None:
            if cmap==None: cmap = get_cmap()
            cmap.set_over(alpha=0.0)
            cmap.set_under(alpha=0.0)
                          
#       Display the image
#       -----------------
        axes([0,0,1,1],frameon=False)
        axis('off')
        im = imshow(X, interpolation=interpolation, cmap=cmap, \
                        aspect='auto', origin='lower', \
                        vmin=vmin, vmax=vmax, **kwopts)
#       All done
#       --------
###        if isinteractive(): show()

#   ..................................................................

    def imshow (self,expr,**kwopts):
        """
        Wrapper around Basemap.imshow() with axis, colorbar and map
        transformations based on the current dimension environment/
        map projection.
        """
        self._imshow (expr,**kwopts)

    def pcolor (self,expr,**kwopts):
        """
        Wrapper around Basemap.pcolor() with axis, colorbar and map
        transformations based on the current dimension environment/
        map projection.
        """
        self._imshow (expr, pcolor=True, **kwopts)

#   ..................................................................

    def _imshow ( self, expr, dh=None, Nx=None, Ny=None, bgim=None,  
                  mpcol=None, bgcol=None, cmap=None, vmin=None, vmax=None, 
                  sub=None, pcolor=False, interpolation='bicubic', 
                  dlat=None, dlon=None, Map=True, cbar_format=None, 
                  **kwopts ):

        """Evaluates GrADS expression and display it as an image
        (or pseudo-color image) on the current basemap."""

#       Check dim environment
#       ---------------------
        if dh==None:
            dh = self.query("dims", Quiet=True)
        if dh.nz>1 or dh.nt>1:
            raise GrADSError, 'Not a horizontal slice; ' + \
                  'expected (nz,nt)=(1,1) but got (%d,%d)'%(dh.nz,dh.nt)

#       Make sure we have a basemap
#       ---------------------------
        if self.map == None or self.transf == False:
            self.basemap(dh=dh)
        m = self.map
        
#       Evaluate GrADS expression
#       -------------------------
        if self.transf: self.cmd('set lon -180 180',Quiet=True)
        Z = self.exp(expr)
        g = Z.grid
        if self.transf: self.cmd('set x %s %s'%dh.x,Quiet=True)

#       Sizes for transformation/blue marble
#       ------------------------------------
        if Ny==None:
            Ny = 600
        if Nx==None:
            aspect = (m.xmax-m.xmin)/(m.ymax-m.ymin)
            Nx = max(450,int(aspect * Ny))

#       Transform data if necessary
#       ---------------------------
        if self.transf:
            if ( abs(g.lon[-1]-180.) < 0.1*(g.lon[1]-g.lon[0]) ):
                g.lon[-1] = 180. # to avoid float point issues
            Z = m.transform_scalar(Z,g.lon,g.lat,Nx,Ny,masked=True)

#       If no background image specified, use Blue Marble if loaded
#       -----------------------------------------------------------
        if (bgim is None) and (self.blue is not None):
            bgim = self.blue_marble(Nx=Nx,Ny=Ny)
            BlueMarble = True
        else:
            BlueMarble = False

#       When vmin/vmax are specified, make the over/under colors
#       transparent
#       -------------------------------------------------------
        if vmin!=None and vmax!=None:
            if cmap==None: cmap = get_cmap()
            cmap.set_over(alpha=0.0)
            cmap.set_under(alpha=0.0)
                          
#       Display the image
#       -----------------
        fig = gcf()
        if sub==None:
            ax = fig.add_axes(self.axis,axisbg=bgcol)
        else:
            ax = fig.add_subplot(sub,axisbg=bgcol)

#       Display background image
#       ------------------------
        if bgim is not None:
            col = 'y'
            if isinstance(bgim,Image.Image):
                m.imshow(_pil_to_array(bgim).astype(float32)/255.) 
            else:
                m.imshow(bgim)  # background image like a Sat Image
        else:
            col = 'k'
        if mpcol == None: mpcol = col

#       Display the image
#       -----------------
        if pcolor:
            lons,lats = m(*meshgrid(g.lon,g.lat))
            im = m.pcolor(lons,lats, Z, 
                          cmap=cmap, vmin=vmin, vmax=vmax,**kwopts )
        else:
            im = m.imshow(Z, interpolation=interpolation, 
                          cmap=cmap, vmin=vmin, vmax=vmax,**kwopts )
        

        # Draw a colorbar
        # ---------------
        self._colorbar(cbar_format,bgcol)

#       Continents
#       ----------
        if Map and (not BlueMarble):
            self.map.drawcoastlines(color=mpcol)
            m.drawmapboundary()

#       Draw lat/lon axis
#       -----------------
        self._labelAxis(dh,dlat=dlat,dlon=dlon)
    
#       All done
#       --------
###     if isinteractive(): show()
        
    def _colorbar(self,cbar_format,bgcol=None):
        """ Draw a colorbar """
        try:    # Jeff made it easier to draw a color bar  
            self.map.colorbar(drawedges=False,format=cbar_format)

        except: # for older vesions of basemap/matplotlib
            draw() # to cope with a bug in MPL > 0.99
            bbox = ax.get_position()
            if type(bbox) is ListType:
                l,b,w,h = bbox        # older mpl < 0.98
            else:
                l,b,w,h = bbox.bounds # mpl >= 0.98
            cax = axes([l+w+0.02, b, 0.04, h],axisbg=bgcol) # setup colorbar axes.
            colorbar(drawedges=False, cax=cax, format=cbar_format) # draw colorbar
            axes(ax)  # make the original axes current again

#   ..................................................................

    def contourf (self,expr,**kwopts):
        """
        Wrapper around Basemap.contourf() with axis, colorbar and map
        transformations based on the current dimension environment/
        map projection.
        """
        self._contourf(expr,**kwopts)

    def contour(self,expr,**kwopts):
        """
        Wrapper around Basemap.contour() with axis, colorbar and map
        transformations based on the current dimension environment/
        map projection.
        """
        self._contourf(expr,cfill=False,clines=True,**kwopts)

#   ..................................................................

    def _contourf ( self, expr, dh=None, bgim=None, 
                    N=None, V=None, cfill=True, clines=False, 
                    mpcol=None, sub=None, dlat=None, dlon=None,
                    Nx=None, Ny=None, Map=True, cbar_format=None, 
                    **kwopts):
        
        """Evaluates GrADS expression and display it as an image
        on the current basemap."""

#       Check dim environment
#       ---------------------
        if dh==None:
            dh = self.query("dims", Quiet=True)
        if dh.nz>1 or dh.nt>1:
            raise GrADSError, 'Not a horizontal slice; ' + \
                  'expected (nz,nt)=(1,1) but got (%d,%d)'%(dh.nz,dh.nt)

#       Make sure we have a basemap
#       ---------------------------
        if self.map == None or self.transf == False:
            self.basemap(dh=dh)
        m = self.map

#       Sizes for transformation/blue marble
#       ------------------------------------
        if Ny==None:
            Ny = 600
        if Nx==None:
            aspect = (m.xmax-m.xmin)/(m.ymax-m.ymin)
            Nx = max(450,int(aspect * Ny))

#       If no background image specified, use Blue Marble if loaded
#       -----------------------------------------------------------
        if (bgim is None) and (self.blue is not None):
            bgim = self.blue_marble(Nx=Nx,Ny=Ny)
            BlueMarble = True
        else:
            BlueMarble = False

#       Evaluate GrADS expression
#       -------------------------
        if self.transf: self.cmd('set lon -180 180',Quiet=True)
        Z = self.exp(expr)
        g = Z.grid
        if self.transf: self.cmd('set x %s %s'%dh.x,Quiet=True)        

#       Setup axis
#       ----------
        fig = gcf()
        if sub==None:
            if cfill:
                ax = fig.add_axes([0.1,0.1,0.75,0.75])
            else:
                ax = fig.add_axes([0.1,0.1,0.8,0.8])
        else:
            ax = fig.add_subplot(sub)

#       Display the contour map
#       -----------------------
        X,Y = m(*meshgrid(g.lon,g.lat))
        if bgim is not None:
            col = 'y'
            if isinstance(bgim,Image.Image):
                m.imshow(_pil_to_array(bgim).astype(float32)/255.) 
            else:
                m.imshow(bgim)  # background image like a Sat Image
        else:
            col = 'k'
        if mpcol == None:  mpcol=col
        if V==None:
            if N==None: N=16
            if clines and cfill:
                cs = m.contour(X,Y,Z,N,linewidths=0.5,colors='k')
            elif clines:
                cs = m.contour(X,Y,Z,N,linewidths=1.25,**kwopts)
                cs.clabel(fmt='%1.0f') 
            if cfill:
                cs = m.contourf(X,Y,Z,N,**kwopts)
        else:
            if clines and cfill:
                cs = m.contour(X,Y,Z,V,linewidths=0.5,colors='k')
            elif clines:
                cs = m.contour(X,Y,Z,V,linewidths=1.25,**kwopts)
                cs.clabel(fmt='%1.0f') 
            if cfill:
                cs = m.contourf(X,Y,Z,V,**kwopts)
            
#       Color Bar
#       ---------
        if cfill:
            self._colorbar(cbar_format)

#       Continents
#       ----------
        if Map and (not BlueMarble):
            m.drawcoastlines(color=mpcol)
            m.drawmapboundary()

#       Draw lat/lon axis
#       -----------------
        self._labelAxis(dh,dlat=dlat, dlon=dlon)
    
#       All done
#       --------
###        if isinteractive(): show()

#   ..................................................................

    def quiver ( self, uexpr, vexpr, dh=None, bgim=None, color=None,
                       mpcol=None, sub=None, dlat=None, dlon=None,
                       Nx=None, Ny=None, vNx=None, vNy=None, Map=True,
                       label=None, labelval=None, labelpos='W', **kwopts):
        
        """Evaluates GrADS expressions with XY vector components
        and display it as on the current basemap."""

#       Check dim environment
#       ---------------------
        if dh==None:
            dh = self.query("dims", Quiet=True)
        if dh.nz>1 or dh.nt>1:
            raise GrADSError, 'Not a horizontal slice; ' + \
                  'expected (nz,nt)=(1,1) but got (%d,%d)'%(dh.nz,dh.nt)

#       Make sure we have a basemap
#       ---------------------------
        if self.map == None or self.transf == False:
            self.basemap(dh=dh)
        m = self.map

#       Sizes for transformation/blue marble
#       ------------------------------------
        if Ny==None:
            Ny = 600
        if Nx==None:
            aspect = (m.xmax-m.xmin)/(m.ymax-m.ymin)
            Nx = max(450,int(aspect * Ny))

#       Sizes for vector transformation
#       -------------------------------
        if vNx==None:
            vNx = dh.nx
        if vNy==None:
            vNy = dh.ny

#       If no background image specified, use Blue Marble if loaded
#       -----------------------------------------------------------
        if (bgim is None) and (self.blue is not None):
            bgim = self.blue_marble(Nx=Nx,Ny=Ny)
            BlueMarble = True
            if color==None: color = 'y'
        else:
            BlueMarble = False

#       Evaluate GrADS expression
#       -------------------------
        if self.transf: self.cmd('set lon -180 180',Quiet=True)
        U = self.exp(uexpr)
        V = self.exp(vexpr)
        g = U.grid
        if self.transf: self.cmd('set x %s %s'%dh.x,Quiet=True)        

#       Setup axis
#       ----------
        fig = gcf()
        if sub==None:
            ax = fig.add_axes([0.1,0.1,0.8,0.8])
        else:
            ax = fig.add_subplot(sub)

        if bgim is not None:
            if mpcol==None: mpcol = 'y'
            m.imshow(bgim)  # background image like a Sat Image

#       Display the vectors
#       -------------------
        if self.transf or (vNx!=dh.nx) or (vNy!=dh.ny):
            U, V, X, Y = m.transform_vector(U, V, g.lon, g.lat, vNx, vNy,
                                            returnxy=True,masked=True)
        else:
            X,Y = m(*meshgrid(g.lon,g.lat))

        self.Q = m.quiver(X,Y,U,V,color=color,**kwopts)
            
#       Color Bar
#       ---------
        if label!=None:
            if labelval == None:
                labelval = sqrt(U*U+V*V).max()
            quiverkey(self.Q,0.15,0.05,labelval,label,
                      coordinates='axes',labelpos=labelpos)
#       Continents
#       ----------
        if Map and (not BlueMarble):
            m.drawcoastlines(color=mpcol)
            m.drawmapboundary()

#       Draw lat/lon axis
#       -----------------
        self._labelAxis(dh,dlat=dlat, dlon=dlon)
    
#       All done
#       --------
###        if isinteractive(): show()

    def imexp (self,imsize=None):
        """Displays GrADS graphics buffer on the Matplotlib canvas.
        On input, size is a tuple with the size of the image, e.g.,
        imsize = (800,600)"""
        imfile = mktemp('.png') 
        if imsize is None:
            self.cmd('gxyat %s'%imfile,Quiet=True)
        else:
            self.cmd('gxyat -x %d -y %d %s'%(imsize[0],imsize[1],imfile),Quiet=True)
        axes([0,0,1,1],frameon=False)
        axis('off')
        imshow(imread(imfile))
        os.remove(imfile)

#   ..................................................................
    def _labelAxis(self,dh,dlat=None,dlon=None):
        """Label axis (private). """

#       Parallels
#       ---------
        if dlat is None:
            LY = dh.lat[1] - dh.lat[0]
            if   LY > 90:  dlat = 30.
            elif LY > 45:  dlat = 15.
            elif LY > 20:  dlat = 5.
            elif LY > 10:  dlat = 2.
            else:          dlat = 1.
        circles = arange(-90.,90.+dlat,dlat).tolist()
        self.map.drawparallels(circles,labels=self.labels)

#       Meridians
#       ---------
        if dlon is None:
            LX = dh.lon[1] - dh.lon[0]
            if   LX > 180:  dlon = 60.
            elif LX >  90:  dlon = 30.
            elif LX >  45:  dlon = 10.
            elif LX >  20:  dlon = 4.
            else:           dlon = 2.
        meridians = arange(-180.,360.+dlon,dlon)
        self.map.drawmeridians(meridians,labels=self.labels)

