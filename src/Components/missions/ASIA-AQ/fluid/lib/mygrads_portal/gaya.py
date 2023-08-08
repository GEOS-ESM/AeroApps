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
This module implements the GaYa class, an extention of the GrADS
client class which adds MayaVi specific methods for 3D visualization.
This class extends GaLab, providing high level methods operating
directly on GrADS expressions (or *fields*) while retaining all the
configurability that MayaVi has to offer.
"""

__version__ = '1.0.0'

from numpy       import mgrid
from galab       import *
from mayavi.mlab import *

class GaYa(GaLab):
    """
    This class extends the GrADS client class GaLab by providing an interface
    to MayaVi. The main design philosophy is to use the
    GrADS dimension environment to derive sensible defaults for the
    high level graphical commands, while retaining most of MayaVi
    configurability.
    """
    
    def __init__ (self, **kwopts):
        """ 
        Contructor. All keyword options (**kwopts) are passed to the
        costructor of the core GrADS client class.
        """
  
        GaLab.__init__(self, **kwopts)

#   ..................................................................

    def surf(self, *args, **kwargs):
        """
        Evaluates GrADS expression and display with the Mayavi
        function *surf*.
        """
        return self._surf(surf, *args, **kwargs)

#   ..................................................................

    def contour_surf(self, *args, **kwargs):
        """
        Evaluates GrADS expression and display with the Mayavi
        function *surf*.
        """
        return self._surf(contour_surf, *args, **kwargs)

#   ..................................................................

    def barchart(self, *args, **kwargs):
        """
        Evaluates GrADS expression and display with the Mayavi
        function *barchart*.
        """
        return self._surf(barchart, *args, **kwargs)

#   ..................................................................

    def _surf(self, what, expr, dh=None, 
              nb_labels=7,xlabel=None,ylabel=None,
              z_axis_visibility=False, # warp_scale='auto',
              **kwopts):
        """
        Evaluates GrADS expression and display with the Mayavi
        function *what*.
        """

        # Check dim environment
        # ---------------------
        if dh==None:
            dh = self.query("dims", Quiet=True)
 #       if dh.nz>1 or dh.nt>1:
 #           raise GrADSError, 'Not a horizontal slice; ' + \
 #                 'expected (nz,nt)=(1,1) but got (%d,%d)'%(dh.nz,dh.nt)

        # Evaluate GrADS expression
        # -------------------------
        Z = self.exp(expr)
        g = Z.grid
        x = g.__dict__[g.dims[1]]
        y = g.__dict__[g.dims[0]]
        X, Y = meshgrid(x,y)
        X, Y, Z = (X.T, Y.T, Z.T)
    
        # Render Surface
        # --------------
        sfc = what(X,Y,Z, #warp_scale=warp_scale,
                   **kwopts)

        # Draw Axes
        # ---------
        if nb_labels>2:
            if xlabel == None: xlabel = g.dims[1]
            if ylabel == None: ylabel = g.dims[0]
            axes(nb_labels=nb_labels,xlabel=xlabel,ylabel=ylabel,
                 z_axis_visibility=z_axis_visibility)

#       All done
#       --------
        return sfc
            
#   ..................................................................

    def contour3d(self, expr, dh=None, 
                  nb_labels=7,xlabel='Longitude',ylabel='Latitude',
                  warp_scale = 1,
                  **kwopts):
        """
        Evaluates GrADS expression and display with the Mayavi
        function *contour3d*.
        """

        # Check dim environment
        # ---------------------
        if dh==None:
            dh = self.query("dims", Quiet=True)
        if dh.nz<2 or dh.nt>1:
            raise GrADSError, 'Not a (X,Y,Z) field; ' + \
                  'expected (nz,nt)=(>1,1) but got (%d,%d)'%(dh.nz,dh.nt)

        # Evaluate GrADS expression
        # -------------------------
        S = self.exp(expr)

        # Assumes an uniform lat/lon grid, use index space for vertical
        # -------------------------------------------------------------
        g = S.grid
        nx, ny, nz = len(g.lon)*1j, len(g.lat)*1j, len(g.lev)
        X, Y, Z = mgrid[g.lon[0]:g.lon[-1]:nx,
                        g.lat[0]:g.lat[-1]:ny,
                        1:warp_scale*nz:nz*1j]
        S = S.T
    
        # Render Surface
        # --------------
        obj = contour3d(X,Y,Z,S,**kwopts)

        # Draw Axes
        # ---------
        if nb_labels>2:
            axes(nb_labels=nb_labels,xlabel=xlabel,ylabel=ylabel)

#       All done
#       --------
        return obj
            
#   ..................................................................

