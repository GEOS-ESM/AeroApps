"""
A Python interface to CPTEC's VegType retrieve Model.
"""

import sys
from numpy      import any

from mxd14      import MxD14_L2
from pyobs.igbp import * 

#--

class IGBP(MxD14_L2):

    def getSimpleVeg(self,Path,nonVeg=GRASSLAND): 
        """
        Given the pathname for the IGBP datasets, attach
    information about the *aggregated* vegetation type:  

         1  Tropical Forests        IGBP  2, 30S < lat < 30N
         2  Extra-Tropical Forests  IGBP  1, 2(lat <=30S or lat >= 30N), 
                                          3, 4, 5
         3  cerrado/woody savanna   IGBP  6 thru  9
         4  Grassland/cropland      IGBP 10 thru 17

     the new attribute name is "veg". Notice that this module also
     defines the following constants:

         TROPICAL = 1
         EXTRA_TROPICAL = 2
         SAVANNA = 3
         GRASSLAND = 4

     corresponding to each aggregated vegetation type."""

#       Vegetation index
#       ----------------
        self.veg = getSimpleVeg(self.lon,self.lat,Path,nonVeg)

#       Assign default flaming fraction, if not set
#       -------------------------------------------
        if self.r_F is None:
            self.r_F = zeros(self.lon.shape)
            b = 1
            for r_F in FLAMING_FRACTION:
                j = (self.veg==b)
                if any(j):
                    self.r_F[j] = r_F
                b += 1
            
#--
    def getDetailedVeg(self,Path): 
        """
        Given the pathname for the IGBP datasets, attach
    information about the vegetation type.
        """
        self.veg = getDetailedVeg(self.lon,self.lat,Path)

if __name__ == "__main__":

    from optparse import OptionParser
    from numpy    import linspace, meshgrid, zeros

    from gfio     import GFIO
    from qfed.mxd14_l3 import __VERSION__ as QFED_VERSION


    # default values for the command line options
    IGBP_DIR = '/nobackup/Emissions/Vegetation/GL_IGBP_INPE'
    QFED_VEG_FILE = 'QFED.vegetation_map.x288_y181.2010.nc4'
    GRID_RESOLUTION = 'c'


    # parse command line options
    parser = OptionParser(usage="Usage: %prog [options]", version='qfed_vegtypes-1.0.b3')

    parser.add_option("-i", "--igbp", 
                      dest="igbp_dir", default=IGBP_DIR,
                      help="directory for IGBP vegetation database (default=%s)" % IGBP_DIR)
    
    parser.add_option("-o", "--output", 
                      dest="qfed_veg", default=QFED_VEG_FILE,
                      help="QFED vegetation map output NetCDF file (default=%s)" % QFED_VEG_FILE)

    parser.add_option("-r", "--resolution", 
                      dest="res", default=GRID_RESOLUTION,
                      help="horizontal resolution: 'a' for 4x5, 'b' for 2x2.5, etc. (default=%s)" % GRID_RESOLUTION)
    
    (options, args) = parser.parse_args()

    if len(args) > 0:
        parser.error("Usage: %prog [options].\nTry '%prog --help' for more information.")
    

    # generate the aggregated vegetation map
    veg_map = getSimpleVegMap(options.igbp_dir, res=options.res)

    f = GFIO()
       
    n_ymd = 20100609
    n_hms = 120000

    v_name  = ['vegetation', ]
    v_title = ['QFED v%3.1f vegetation map' % QFED_VERSION, ]
    v_units = ['', ]

    title = 'QFED v%3.1f Vegetation Map' % QFED_VERSION
    source = 'NASA/GSFC/GMAO GEOS-5 Aerosol Group'
    contact = ('%s; %s') % ('arlindo.dasilva@nasa.gov', 'anton.darmenov@nasa.gov')

    # create the output file
    try:
        f.create(options.qfed_veg,
                 v_name,
                 n_ymd,
                 n_hms,
                 lon=veg_map['lon'],
                 lat=veg_map['lat'],
                 vtitle=v_title,
                 vunits=v_units,
                 title=title,
                 source=source,
                 contact=contact)
    except:
        print "ERROR: Could not create QFED vegetation file."
        sys.exit(1)

    # write the data
    try:
        f.write(v_name[0], n_ymd, n_hms, veg_map['data'])
    except:
        print "ERROR: Could not save the QFED vegetation file."

        try:
            f.close()
        except:
            pass

        sys.exit(2)

    # close the file
    try: 
        f.close()
    except: 
        pass

    print "INFO: The QFED v%3.1f vegetation map was saved to '%s'\n" % (QFED_VERSION, options.qfed_veg)

