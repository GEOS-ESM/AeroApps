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
    n_split = 1

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

    parser.add_option("-y", "--Ny", dest="Ny", default=None,
              type='int',
              help="Number of NORTH-SOUTH pixels (default:same as Nx)" )

    parser.add_option("-s", "--n_split", dest="n_split", default=n_split,
              type='int',
              help="Number of time slices for each 30 min global image (default:%d)"%n_split )

    parser.add_option("-p", "--prime",
                      action="store_true", dest="prime",
                      help="Input image starts at 17.5W; default is date line")

    parser.add_option("-I", "--interpolate",
                      action="store_true", dest="interp",
                      help="Time interpolate in between images")

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

    if options.Ny is None:
        options.Ny = options.Nx

    if len(Images)==1 or options.n_split==1:
        options.interp = False # no time interpolation in this case.

    # Time steps
    # ----------
    DT = timedelta(minutes=30) # global images are 30 min (hardwired)
    dt = DT / options.n_split
        
    # Upstream image
    # --------------
    if options.interp:
        M = len(Images)-1
        RGB2 = imread(Images[0])
    else:
        M = len(Images)

    # Loop over images
    # ----------------
    for i in range(len(Images)):

        imFile = Images[i]

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

        # Time interpolating
        # ------------------
        if options.interp:
            RGB1 = RGB2
            RGB2 = imread(Images[i+1])
        else:
            RGB1 = imread(imFile)

        nlats, nlons, nc = RGB1.shape

        # Define lat/lon grid that image spans (projection='cyl').
        # --------------------------------------------------------
        dlon = 360./float(nlons)
        dlat = 180./float(nlats-1)
        lons = arange(-180,180.,dlon)
        lats = arange(-90.,90+dlat,dlat)
        I = range(nlons)
        if options.prime: # notice 17.5W, not quite prime meridian
            Lons = arange(-17.5,360.-17.5,dlon)
            J = Lons>180
            Lons[J] = Lons[J] - 360.
            I = Lons.argsort() # this will do a lon swap so that array goes from -180 to 180.
       

        # Get image time from file name
        # -----------------------------
        tokens = os.path.basename(imFile).split('_')
        nymd = tokens[5]
        hhmm = tokens[6]
        Y, M, D = int(nymd[0:4]), int(nymd[4:6]), int(nymd[6:])
        h, m = int(hhmm[0:2]),int(hhmm[2:4])
        tyme0 = datetime(Y,M,D,h,m)

        # Time slices (if so desired for smoother animation)
        # --------------------------------------------------
        for n in range(options.n_split):

            # This image time
            # ---------------
            tyme = [tyme0 + n * dt,]
            
            t = tyme[0]
            Y,M,D,h,m = (t.year,t.month,t.day,t.hour,t.minute)
            outfile = outFile[:-18] + '%4d%02d%02d_%02d%02dz.png'%(Y,M,D,h,m) 

            if os.path.exists(outfile):
                print "[x] Output file exists, skipping %s"%outfile
                continue

            if options.verbose:
                print "    Working on %s"%outfile

            # Get central (lon,lat) from time
            # -------------------------------
            lon_0, lat_0 = getCenterPoint(tyme)
            
            # Interpolate global image in time
            # --------------------------------
            if options.interp:
                a = float(n) / float(options.n_split)
                RGB = (1-a) * RGB1 + a * RGB2
            else:
                RGB = RGB1

            # Instantiate basemap transform for this time
            # Not sure why, but we need lat_0 = - lat_0
            # -------------------------------------------
            m =  Basemap(projection='ortho',resolution=None,
                         lon_0 = lon_0[0], lat_0 = -lat_0[0],
                         rsphere=(rsphere[0]+rsphere[1])/2)

            # Interpolate from global latlon to ortho disk
            # --------------------------------------------
            rgb = ma.zeros((options.Ny,options.Nx,nc))
            for k in range(nc):
                if options.verbose:
                    print "    - Remapping Layer ", k, options.Nx, options.Ny, lon_0[0], lat_0[0]
                rgb[:,:,k] = m.transform_scalar(RGB[:,I,k],lons,lats,
                                                options.Nx, options.Ny,
                                                masked=True)

            rgb[rgb.mask] = 0

            # Save the image
            # --------------
            imsave(outfile,rgb)
        
