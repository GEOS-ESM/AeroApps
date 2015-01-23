#!/usr/bin/env python

"""
    Utility to project a GEOPS-5 global image file on a geostationary
    sector.
"""

import os

from netCDF4 import Dataset

from numpy import zeros, ones, arange, array, ma, linspace, meshgrid
from mpl_toolkits.basemap import Basemap, interp

from matplotlib.image import imread, imsave

from optparse        import OptionParser
from datetime        import datetime, timedelta
from dateutil.parser import parse as isoparser

from inclination import getCenterPoint

# Earth model
# -----------
rsphere=(6378137.00,6356752.3142)
satellite_height = 35785831.0

# -------------------------- M A I N -----------------------------
if __name__ == "__main__":
    
    
    # DSCOVR based defaults
    # ---------------------
    Nx = 2880 # number of X pixels
    outDir = '.'
    rPat='_globe_,_dscovr_'

#   Parse command line options
#   --------------------------
    parser = OptionParser(usage="Usage: %prog [OPTIONS] geos5_global_image(s)",
                          version='1.0.0' )

    parser.add_option("-d", "--outDir", dest="outDir", default=outDir,
              help="Output directory name (default=%s)"\
                          %outDir )

    parser.add_option("-r", "--rPat", dest="rPat", default=rPat,
              help="Output filename replacement pattern (default=%s)"\
                          %rPat )

    parser.add_option("-x", "--Nx", dest="Nx", default=Nx,
              type='int',
              help="Number of EAST-WEST pixels (default=%d)"\
                          %Nx )

    parser.add_option("-y", "--Ny", dest="Ny", default=Nx,
              type='int',
              help="Number of NORTH-SOUTH pixels (default is same as Nx)" )

    parser.add_option("-f", "--force",
                      action="store_true", dest="force",
                      help="Overwrite existing output files")

    parser.add_option("-v", "--verbose",
                      action="store_true", dest="verbose",
                      help="Verbose mode")

    # Parse arguments
    # ---------------
    (options, args) = parser.parse_args()
    if len(args) < 1:
        parser.error("must have 1 argument: geos5_global_image_file")

    Images = sorted(args)
    p1, p2 = options.rPat.split(',')

    # Loop over images
    # ----------------
    for imFile in Images:

        # Get image time from file name
        # -----------------------------
        tokens = imFile.split('_')
        nymd = tokens[5]
        hhmm = tokens[6]
        Y, M, D = int(nymd[0:4]), int(nymd[4:6]), int(nymd[7:])
        h, m = int(hhmm[0:2]),int(hhmm[2:4])
        tyme = [ datetime(Y,M,D,h,m), ]
        
        # Get centrol (lon,lat) from time
        # -------------------------------
        lon_0, lat_0 = getCenterPoint(tyme)
            
        # Instantiate basemap transform for this time
        # -------------------------------------------
        m =  Basemap(projection='ortho',resolution=None,
                     lon_0 = lon_0, lat_0 = lat_0,
                     rsphere=(rsphere[0]+rsphere[1])/2)

        # Form output file
        # ----------------
        outFile = options.outDir + '/' + os.path.basename(imFile).replace(p1,p2)

        # Unless force, do not overwrite files
        # ------------------------------------
        if not options.force:
            if os.path.exists(outFile):
                print "[x] Output file exists, skipping %s"%outFile
                continue

        if options.verbose:
            print "[ ] Reading %s"%os.path.basename(imFile)

        # Use PIL to read image
        # ---------------------
        RGB = imread(imFile)
        nlats, nlons, nc = RGB.shape

        # Define lat/lon grid that image spans (projection='cyl').
        # --------------------------------------------------------
        dlon = 360./float(nlons)
        dlat = 180./float(nlats-1)
        lons = arange(-180,180.,dlon)
        lats = arange(-90.,90+dlat,dlat)

        # Interpolate from global latlon to ortho disk
        # --------------------------------------------
        rgb = ma.zeros((options.Ny,options.Nx,nc))
        for k in range(nc):
            if options.verbose:
                print "   - Remapping Layer ", k
            rgb[:,:,k] = m.transform_scalar(RGB[:,:,k],lons,lats,
                                            options.Nx, options.Ny,
                                            masked=True)
        # Save the image
        # --------------
        if options.verbose:
            print "   - Writing %s"%outFile
        imsave(outFile,rgb)
        
