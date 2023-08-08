#--------------------------------------------------------------------------
#
#    Copyright (C) 2006-2008 by Arlindo da Silva <dasilva@opengrads.org>
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
Define additional color tables.
"""

__version__ = '1.2.0'

from numpy import log10, exp, log, arange, ma
import numpy as npy
from matplotlib        import rcParams, cbook, cm
from matplotlib.colors import Colormap, makeMappingArray
from types import *

LUTSIZE = rcParams['image.lut']

class ColormapAlpha(Colormap):

    def _call_(self, X, alpha=1, bytes=None):
        """
        OBSOLETE: to be removed
        X is either a scalar or an array (of any dimension).
        If scalar, a tuple of rgba values is returned, otherwise
        an array with the new shape = oldshape+(4,). If the X-values
        are integers, then they are used as indices into the array.
        If they are floating point, then they must be in the
        interval (0.0, 1.0).
          *alpha* must be a scalar.
          *bytes* not used yet.
        """

        if not self._isinit: self._init()
        
        if alpha==1 or alpha==None:
            alpha = self._alpha[:]
        else:
            alpha = min(alpha, 1.0) # alpha must be between 0 and 1
            alpha = max(alpha, 0.0)
            alpha = alpha * self._alpha

        self._lut[:-3,3] = alpha[:]

        mask_bad = None
        if not cbook.iterable(X):
            vtype = 'scalar'
            xa = npy.array([X])
        else:
            vtype = 'array'
            xma = ma.asarray(X)
            xa = xma.filled(0)
            mask_bad = ma.getmask(xma)
        if xa.dtype.char in npy.typecodes['Float']:
            npy.putmask(xa, xa==1.0, 0.9999999) #Treat 1.0 as slightly less than 1.
            xa = (xa * self.N).astype(int)
        # Set the over-range indices before the under-range;
        # otherwise the under-range values get converted to over-range.
        npy.putmask(xa, xa>self.N-1, self._i_over)
        npy.putmask(xa, xa<0, self._i_under)
        if mask_bad is not None and mask_bad.shape == xa.shape:
            npy.putmask(xa, mask_bad, self._i_bad)
        rgba = self._lut[xa]
        if vtype == 'scalar':
            rgba = tuple(rgba[0,:])
        return rgba


class SegmentedColormap(ColormapAlpha):
    """Colormap objects based on lookup tables using generic segments.

    The lookup transfer function is a simple linear function between
    defined intensities. There is no limit to the number of segments
    that may be defined. Though as the segment intervals start containing
    fewer and fewer array locations, there will be inevitable quantization
    errors.

    This customized version allows one to define the alpha channel as well.

    """
    def __init__(self, name, segmentdata, N=256, scale=None, reverse=False):
        """Create color map from linear mapping segments
        segmentdata argument is a dictionary with a red, green and blue
        entries. Each entry should be a list of x, y0, y1 tuples.
        See makeMappingArray for details. Optionally, the x coordinate
        in [0,1] is mapped thru scale(x).
        """
        self.monochrome = False  # True only if all colors in map are identical;
                                 # needed for contouring.
        ColormapAlpha.__init__(self, name, N)
        self._reverse = reverse
        self._scale = scale
        self._segmentdata = segmentdata
            
    def _init(self):
        """Initialize from a tradition segment data or from an
           existing colormap redefining the intrinsic alpha channel.
        """

#       Inherit segment data from an existing colormap,
#       possibly redefining alpha channel
        if self._segmentdata.has_key('cmap'):
            if self._segmentdata.has_key('alpha'):
                alphad = self._segmentdata['alpha'] # save this
            else:
                alphad = None
            cmap = self._segmentdata['cmap']
            if not cmap._isinit: cmap._init()
            self._segmentdata = cmap._segmentdata.copy()
            if alphad!=None:
                self._segmentdata['alpha'] = alphad # revised alpha

#       Scale the segment data if requested
#       -----------------------------------
        if self._scale!=None:
            scale = lambda x: min(max(self._scale(x),0.0),1.0)
            data_s = {}
            for key, val in self._segmentdata.iteritems():
                if not isinstance(val, list): continue 
                valnew = [ (scale(a), b, c) for a, b, c in val ]
                print valnew
                data_s[key] = valnew
            self._segmentdata = data_s

#       Reverse the segment data if requested
#       -------------------------------------
        if self._reverse:
            data_r = {}
            for key, val in self._segmentdata.iteritems():
                valnew = [(1.0-a, b, c) for a, b, c in list(reversed(val))]
                data_r[key] = valnew
            self._segmentdata = data_r

#       Go for the LUT
#       --------------
        self._lut = npy.ones((self.N + 3, 4), npy.float)
        self._lut[:-3, 0] = makeMappingArray(self.N, self._segmentdata['red'])
        self._lut[:-3, 1] = makeMappingArray(self.N, self._segmentdata['green'])
        self._lut[:-3, 2] = makeMappingArray(self.N, self._segmentdata['blue'])

#       RGB is always there by alpha may not, so we check for it
#       --------------------------------------------------------
        if self._segmentdata.has_key('alpha'):
            self._lut[:-3, 3] = makeMappingArray(self.N, self._segmentdata['alpha'])

#       Save intrinsic alpha
#       --------------------
        self._alpha = self._lut[:-3, 3][:] 

        self._isinit = True
        self._set_extremes()

#............................................................................

def nlog(x,a):
    if a<=0:
        raise ValueError, 'Expected a>0 but got a=%f'%a
    return log(a*x+1.0)/log(a+1.0)
def nexp(x,a):
    if a<=0:
        raise ValueError, 'Expected a>0 but got a=%f'%a
    return (exp(x * log(a+1.0)) - 1.0) / a

#
# Define additional transparent color tables.
#

datad = {}

# From Jeff DLB: Dust AOT - use (vmax,vmin) = (0,2)   
# red = [0.02, 0.2], green = [0.2, 2.0], blue = [1.6, 2.0], alpha = [0.02, 1.0]
# red = [0.01, 0.1], green = [0.1, 1.0], blue = [0.8, 1.0], alpha = [0.02, 1.0]
# tickname = ['0.02', '0.05', '0.1', '0.2', '0.5', '1.0', '2.0']

datad['dust'] = {   'red': ((0.00, 0.0, 0.0),
                            (0.01, 0.0, 0.0),
                            (0.10, 1.0, 1.0),
                            (1.00, 1.0, 1.0)),
                  'green': ((0.00, 0.0, 0.0),
                            (0.10, 0.0, 0.0),
                            (1.00, 1.0, 1.0)),
                   'blue': ((0.00, 0.0, 0.0),
                            (0.80, 0.0, 0.0),
                            (1.00, 1.0, 1.0)),
                  'alpha': ((0.00, 0.0, 0.0),
                            (0.01, 0.0, 0.0),
                            (0.50, 1.0, 1.0),
                            (1.00, 1.0, 1.0)) 
                }

datad['aot_du'] = datad['dust']

# Seasalt
#   red = [1,1]
#   green = [1,1]
#   blue = [0.01, 1.0]
#   alpha = [0.01, 0.1]
#   tickname = ['0.01', '0.02', '0.05', '0.1', '0.2', '0.5', '1.0'] ;;LOG

datad['aot_ss'] = {'red': ((0.00, 0.0, 0.0),
                            (1.00, 0.0, 0.0)),
                  'green': ((0.00, 0.0, 0.0),
                           (1.00, 0.0, 0.0)),
                   'blue': ((0.00, 0.0, 0.0),
                            (0.01, 0.0, 0.0),
                            (1.00, 1.0, 1.0)),
                  'alpha': ((0.00, 0.0, 0.0),
                            (0.01, 0.0, 0.0),
                            (0.50, 1.0, 1.0),
                            (1.00, 1.0, 1.0)) 
                }

# Sulfates
# red = [0.1, 1.0]
# green = [0.01, 1.0]
# blue = [0.8, 1.0]
# alpha = [0.01, 0.2]
# tickname = ['0.01', '0.02', '0.05', '0.1', '0.2', '0.5', '1.0'] 

datad['aot_su'] = {'red': ((0.00, 0.0, 0.0),
                            (0.10, 0.0, 0.0),
                            (1.00, 1.0, 1.0)),
                  'green': ((0.00, 0.0, 0.0),
                            (0.01, 0.0, 0.0),
                            (1.00, 1.0, 1.0)),
                   'blue': ((0.00, 0.0, 0.0),
                            (0.80, 0.0, 0.0),
                            (1.00, 1.0, 1.0)),
                  'alpha': ((0.00, 0.0, 0.0),
                            (0.01, 0.0, 0.0),
                            (0.50, 1.0, 1.0),
                            (1.00, 1.0, 1.0)) 
                }

# Organic Carbon AOT
#   red = [0.01, 0.1]
#   green = [0.1, 1.0]
#   blue = [0.01, 1.0]
#   alpha = [0.01, 0.1]
#   tickname = ['0.01', '0.02', '0.05', '0.1', '0.2', '0.5', '1.0'] 

datad['aot_oc'] = {'red': ((0.00, 0.0, 0.0),
                            (0.01, 0.0, 0.0),
                            (0.10, 1.0, 1.0),
                            (1.00, 1.0, 1.0)),
                  'green': ((0.00, 0.0, 0.0),
                            (0.10, 0.0, 0.0),
                            (1.00, 1.0, 1.0)),
                   'blue': ((0.00, 0.0, 0.0),
                            (0.01, 0.0, 0.0),
                            (1.00, 1.0, 1.0)),
                  'alpha': ((0.00, 0.0, 0.0),
                            (0.01, 0.0, 0.0),
                            (0.50, 1.0, 1.0),
                            (1.00, 1.0, 1.0)) 
                }

# Black Carbon
#   red = [0.002, 0.2]  --> 0.01 1.
#   green = [0.02, 0.2] --> 0.1  1.
#   blue = [0.002, 0.02] --> 0.01 0.1
#   alpha = [0.002, 0.02]
#   tickname = ['0.002', '0.005', '0.01', '0.02', '0.05', '0.1', '0.2']

datad['aot_bc'] = {'red': (( 0.,    0., 0.),
                            (0.01,  0., 0.),
                            (1.,    1., 1.)
                            ),
                  'green': ((0.,    0., 0.),
                            (0.1,   0., 0.),
                            (1.,    1., 1.)
                            ),
                   'blue': ((0.,    0., 0.),
                            (0.01,  0., 0.),
                            (0.1,   1., 1.),
                            (1.,    1., 1.)
                            ),
                 'alpha': ((0.00, 0.0, 0.0),
                            (0.001, 0.0, 0.0),
                            (0.50, 1.0, 1.0),
                            (1.00, 1.0, 1.0)
                            ) 
                }

alphad = ((0.00, 0.0, 0.0),
          (0.01, 0.0, 0.0),
          (0.50, 1.0, 1.0),
          (1.00, 1.0, 1.0))


# Customizing standard color maps, opaque

datao = dict()
datao['Hot'] = { 'cmap': cm.hot, }
datao['Jet'] = { 'cmap': cm.jet, }
datao['Blues'] = { 'cmap': cm.Blues, }
datao['Oranges'] = { 'cmap': cm.Oranges, }
datao['Pinks'] = { 'cmap': cm.pink, }
datao['Reds'] = { 'cmap': cm.Reds, }
datao['Greens'] = { 'cmap': cm.Greens, }
datao['Greys'] = { 'cmap': cm.Greys, }
datao['Copper'] = { 'cmap': cm.copper, }
datao['Bone'] = { 'cmap': cm.bone, }
datao['Purples'] = { 'cmap': cm.Purples, }
datao['Spectral'] = { 'cmap': cm.spectral, }
datao['Earth'] = { 'cmap': cm.gist_earth, }
datao['Stern'] = { 'cmap': cm.gist_stern, }
datao['Flag'] = { 'cmap': cm.flag, }
datao['Terrain'] = { 'cmap': cm.terrain, }
datao['Seismic'] = { 'cmap': cm.seismic, }


# Customizing standard color maps, transparent

datad['hot'] = { 'cmap': cm.hot, 'alpha': alphad }
datad['jet'] = { 'cmap': cm.jet, 'alpha': alphad }
datad['blues'] = { 'cmap': cm.Blues, 'alpha': alphad }
datad['oranges'] = { 'cmap': cm.Oranges, 'alpha': alphad }
datad['pinks'] = { 'cmap': cm.pink, 'alpha': alphad }
datad['reds'] = { 'cmap': cm.Reds, 'alpha': alphad }
datad['greens'] = { 'cmap': cm.Greens, 'alpha': alphad }
datad['greys'] = { 'cmap': cm.Greys, 'alpha': alphad }
datad['copper'] = { 'cmap': cm.copper, 'alpha': alphad }
datad['bone'] = { 'cmap': cm.bone, 'alpha': alphad }
datad['purples'] = { 'cmap': cm.Purples, 'alpha': alphad }
datad['spectral'] = { 'cmap': cm.Spectral, 'alpha': alphad }
datad['earth'] = { 'cmap': cm.gist_earth, 'alpha': alphad }
datad['stern'] = { 'cmap': cm.gist_stern, 'alpha': alphad }
datad['flag'] = { 'cmap': cm.flag, 'alpha': alphad }
datad['terrain'] = { 'cmap': cm.terrain, 'alpha': alphad }
datad['seismic'] = { 'cmap': cm.seismic, 'alpha': alphad }

datad['hot_r'] = { 'cmap': cm.hot_r, 'alpha': alphad }
datad['jet_r'] = { 'cmap': cm.jet_r, 'alpha': alphad }
datad['blues_r'] = { 'cmap': cm.Blues_r, 'alpha': alphad }
datad['oranges_r'] = { 'cmap': cm.Oranges_r, 'alpha': alphad }
datad['pink_r'] = { 'cmap': cm.pink_r, 'alpha': alphad }
datad['reds_r'] = { 'cmap': cm.Reds_r, 'alpha': alphad }
datad['greens_r'] = { 'cmap': cm.Greens_r, 'alpha': alphad }
datad['greys_r'] = { 'cmap': cm.Greys_r, 'alpha': alphad }
datad['copper_r'] = { 'cmap': cm.copper_r, 'alpha': alphad }
datad['bone_r'] = { 'cmap': cm.bone_r, 'alpha': alphad }
datad['purples_r'] = { 'cmap': cm.Purples_r, 'alpha': alphad }
datad['spectral_r'] = { 'cmap': cm.Spectral_r, 'alpha': alphad }
datad['earth_r'] = { 'cmap': cm.gist_earth_r, 'alpha': alphad }
datad['stern_r'] = { 'cmap': cm.gist_stern_r, 'alpha': alphad }
datad['flag_r'] = { 'cmap': cm.flag_r, 'alpha': alphad }
datad['terrain_r'] = { 'cmap': cm.terrain_r, 'alpha': alphad }
datad['seismic_r'] = { 'cmap': cm.seismic_r, 'alpha': alphad }

#datad[''] = { 'cmap': cm., 'alpha': alphad }

# Define direct, reversed and scaled color tables (transparent)
log_scale = lambda x: nexp(x,10.0) 
for name, data in datad.iteritems():
    locals()[name] = SegmentedColormap(name, data, LUTSIZE)
#    name_r = name+'_r'
#    locals()[name_r] = SegmentedColormap(name_r, data, LUTSIZE, reverse=True )
    name_l = name+'_l'
    locals()[name_l] = SegmentedColormap(name_l, data, LUTSIZE, scale=log_scale)

# Define direct, reversed and scaled color tables (opaque)
for name, data in datao.iteritems():
    #locals()[name] = SegmentedColormap(name, data, LUTSIZE)
    name_r = name+'_r'
    locals()[name_r] = SegmentedColormap(name_r, data, LUTSIZE, reverse=True )
    name_l = name+'_l'
    locals()[name_l] = SegmentedColormap(name_l, data, LUTSIZE, scale=log_scale)
    name_rl = name+'_r_l'
    locals()[name_rl] = SegmentedColormap(name_rl, data, LUTSIZE, reverse=True,scale=log_scale)


"""
From JDLB IDL sources:

Dust AOT
--------

   ;; use log10 scaling
   scaling = 'log'
   ;; Range: 0.01 - 1.0
   red = [0.01, 0.1]
   green = [0.1, 1.0]
   blue = [0.8, 1.0]
   alpha = [0.01, 0.5]
   tickname = ['0.01', '0.02', '0.05', '0.1', '0.2', '0.5', '1.0'] ;;LOG
   ;; Range: 0.02 - 2.0
   red = [0.02, 0.2]
   green = [0.2, 2.0]
   blue = [1.6, 2.0]
   alpha = [0.02, 1.0]
   tickname = ['0.02', '0.05', '0.1', '0.2', '0.5', '1.0', '2.0'] ;;LOG

Black Carbon AOT
----------------
 'BCEXTTAU': begin
   title = 'Black Carbon Aerosol Optical Thickness'
   multiplier = 1D
   scaling = 'log'
   units = ''
   style = 'smooth'
   blue = [0.002, 0.02]
   green = [0.02, 0.2]
   red = [0.002, 0.2]
   alpha = [0.002, 0.02]
   tickname = ['0.002', '0.005', '0.01', '0.02', '0.05', '0.1', '0.2']


Organic Carbon AOT
-----------------
 'OCEXTTAU': begin
   title = 'Organic Carbon Aerosol Optical Thickness'
   red = [0.01, 0.1]
   green = [0.1, 1.0]
   blue = [0.01, 1.0]
   alpha = [0.01, 0.1]
   tickname = ['0.01', '0.02', '0.05', '0.1', '0.2', '0.5', '1.0'] ;;LOG

Organics AOT
------------
 'CCEXTTAU': begin
   title = 'Organic+Black Carbon AOT'
   red = [0.01, 0.1]
   green = [0.1, 1.0]
   blue = [0.01, 1.0]
   alpha = [0.01, 0.1]
   tickname = ['0.01', '0.02', '0.05', '0.1', '0.2', '0.5', '1.0'] ;;LOG

Total AOT
---------
 'TTEXTTAU': begin
   title = 'Total AOT (RGB=dust,sulfates+carbon,salt)'
   multiplier = 1D ;;
   units = ''
   scaling = 'log'
   style = 'rgba' ;;NOTE: this dataset handled differently--see crave.pro
   red = [0.01, 1.0]
   green = [0.01, 1.0]
   blue = [0.01, 1.0]
   alpha = [0.01, 0.1]
   tickname = ['0.01', '0.02', '0.05', '0.1', '0.2', '0.5', '1.0'] ;;LOG

Sea-salt
--------
 'SS-TOTAL': begin ;;linear SS for Total AOT RGB portrayal
   title = 'Sea Salt Total AOT'
   red = [1,1]
   green = [1,1]
   blue = [0.01, 1.0]
   alpha = [0.01, 0.1]
   tickname = ['0.01', '0.02', '0.05', '0.1', '0.2', '0.5', '1.0'] ;;LOG
 end

Sulfate AOT
-----------
 'SUEXTTAU': begin
   title = 'Sulfate Aerosol Optical Thickness'
   red = [0.1, 1.0]
   green = [0.01, 1.0]
   blue = [0.8, 1.0]
   alpha = [0.01, 0.2]
   tickname = ['0.01', '0.02', '0.05', '0.1', '0.2', '0.5', '1.0'] ;;LOG


CO2 - surface
-------------
 'CO2SC001': begin
   title = 'CO2 Surface Concentration'
   red = [360, 390]
   green = [0, 0]
   blue = [360, 390]
   ;; min and max values of transparency ramp
   alpha = [360, 370]
   tickname = ['360', '370', '380', '390']

CO Column
---------
 'COCL001': begin
   title = 'Carbon Monoxide - Global'
   legendtitle = 'CO Column Density'
   multiplier = 1D/0.03 / 1e4 * 6.02e23
   units = 'molec/cm!A2!N'
   scaling = 'log'
   style = 'smooth'
   ;; log scaling
   scaling = 'log'
   red = [1.0e17, 1.0e19, 1.0e19, 0, 255, 255]
   green = [1.0e15, 1.0e16, 1.0e17, 0, 255, 0]
   blue = [1.0e16, 1.0e17, 1.0e18, 0, 255, 0]
   alpha = [1.0e18, 2.0e18, 2.0e18, 0, 255, 255]
   tickname = ['1e15', '1e16', '1e17', '1e18', '1e19']
 end

CO - surface
------------
 'COSC001': begin
   title = 'CO Surface Concentration'
   red = [100.0, 1000.0]
   green = [800.0, 1000.0]
   blue = [50.0, 500.0]
   alpha = [50.0,500.0]
   tickname = ['100', '200', '500', '1000']
 end

Ozone column
------------
 'O3DU': begin
   title = 'Ozone Total Column'
   ;; min, max values to display in each color
   red = [1000, 1000]
   green = [200, 500]
   blue = [200, 500]
   ;; min and max values of transparency ramp
   alpha = [200, 400]
   ;;TO DO: the following 6-parameter portrayal looks great but is not yet
   ;;supported by getlegend.pro.
   ;;alpha = [200, 300, 300, 0, 230, 230]
   tickname = ['200', '300', '400', '500']
   ;; The preceding portrayal is too faint--try to get more detail
   style = 'IDL22'
   red = [200, 400]
   green = red ;;ignored for IDL CT style
   blue = blue ;;ignored for IDL CT style
   alpha = [200, 350]
   tickname = ['200', '300', '400']

"""
