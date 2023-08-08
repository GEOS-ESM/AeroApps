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
This package implements a Python_ interface to GrADS_ by means of
bi-directional pipes. The following modules are provided:

gacore 
  The module *gacore* provides the basic GrADS client class which
  allows you to start GrADS, send commands to it and to retrieve the
  text output produced by the *grads* application in response to such
  command. This is a Pure Python module, although it requires the
  GrADS binaries to have been installed on your system. This module
  defines the class *GaCore*.

ganum 
  If you have NumPy installed, the module *ganum* will be loaded. This
  module defines class *GaNum* which extends the GrADS client class
  *GaCore* by providing methods for exchanging n-dimensional NumPy
  array data between Python and GrADS. It also provides methods for
  computing EOF and least square estimation.

galab 
  If PyLab/Matplotlib/Basemap is available, the module *galab* is
  loaded.  This module defines class *GaLab* which etends class
  *GaNum* with Matplotlib/Basemap specific methods for contours,
  images and other graphical functionality. This class provides high
  level methods operating directly on GrADS expressions (or *fields*)
  while retaining all the configurability that Matplotlib has to
  offer.

gaya
  If MayaVi is available, the module *gaya* is loaded. This module
  defines class *GaYa* with MayaVi specific methods for 3D
  visualization.  This class provides high level methods operating
  directly on GrADS expressions (or *fields*) while retaining all the
  configurability that MayaVi has to offer.

gahandle
  This module provides a simple container class to collect output for
  query() operations.

gacm
  This modules provides additional colormaps, as well as an extension
  of the *Colormaps* class which allows for the definition of color
  look-up takes with an alpha channel. It also extends the
  *LinearSegmentedColormap* with the ability of create derived color
  tables which are either reversed or provide an arbitrary scaling by
  means of a lambda function.

numtypes 
  This module defines GaField, a class for representing GrADS
  variables in Python. It consists of a NumPy masked array with a
  *grid* containing coordinate/dimension information attached to it.

The generic class *GrADS* is an alias to *GaLab* when module *galab*
is sucessfully loaded (usually when NumPy, Matplotlib and PIL are
available.) Otherwise, *GrADS* is an alias to *GaNum* if module
*ganum* can be sucessfully loaded (usually when NumPy is
available). Failing that, the class *GrADS* becomes an alias to
*GaCore* which only relies on the Python standard library. The
following boolean attributes can be used to determine hich
functionality is available: HAS_GALAB, HAS_GANUM, HAS_GACM.

"""

__version__ = '1.2.0'

#
# Import the most capable moodule we have the dependencies for.
# Notice that galab imports all of ganum which in turn imports
# all of gacore.
#

HAS_GACORE = False
HAS_GANUM  = False
HAS_GALAB  = False
HAS_GAYA   = False

try:
    import galab
    GrADS  = galab.GaLab
    GaCore = galab.GaCore
    GaLab  = galab.GaLab
    GaNum  = galab.GaNum
    HAS_GALAB  = True
    HAS_GANUM  = True
    HAS_GACORE = True
    try:
        import gaya
        GrADS = gaya.GaYa
        GaYa  = gaya.GaYa
        HAS_GAYA = True
    except:
        HAS_GAYA = False
except:
    try:
        import ganum
        GrADS  = ganum.GaNum
        GaCore = ganum.GaCore
        GaNum  = ganum.GaNum
        HAS_GANUM  = True
        HAS_GACORE = True
    except:
        import gacore
        GrADS  = gacore.GaCore
        GaCore = gacore.GaCore
        HAS_GACORE = True
try:
    import gacm
    HAS_GACM = True
except:
    HAS_GACM = False

from gacore   import GrADSError
if ( HAS_GANUM ) :
    from numtypes import GaGrid, GaField
