#!/usr/bin/env python3

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
from dateutil.parser import parse         as isoparser

# Earth model
# -----------
rsphere=(6378137.00,6356752.3142)
satellite_height = 35785831.0

class MyBasemap(Basemap):

    def set_bbox(self,bbox):
        """
        Set bounding box for transform operations.
        """
        self.bbox = bbox

    def makegrid(self,nx,ny,returnxy=False,setCache=False):
        """
        return arrays of shape (ny,nx) containing lon,lat coordinates of
        an equally spaced native projection grid.

        If ``returnxy = True``, the x,y values of the grid are returned also.

        *** This fixes cropping bug in the original makegrid. ***

        """

        # For efficiency, cache coordinates
        # ---------------------------------
        try:
            cache = self.makegrid_cache 
        except:
            self.makegrid_cache = False 
            cache = False

        if cache:
            lons, lats, x, y = self.makegrid_grid 
            if returnxy:
                return (lons,lats,x,y)
            else:
                return (lons,lats)
 
        # x and y coordinates
        # -------------------
        x = linspace(self.bbox[0],self.bbox[2],nx,endpoint=True)
        y = linspace(self.bbox[1],self.bbox[3],ny,endpoint=True)
        y = y[-1::-1]
        X, Y = meshgrid(x,y)
        lons,      lats      = (ma.zeros(X.shape), ma.zeros(X.shape))
        lons[:,:], lats[:,:] = self(X,Y,inverse=True)
        lons.set_fill_value(1e30)
        lats.set_fill_value(1e30)
        lons.mask, lats.mask = (lons>1e29, lons>1e29)
 
        if setCache:
            self.makegrid_cache = True
            self.makegrid_grid = (lons,lats,x,y)
            
        if returnxy:
            return (lons,lats,x,y)
        else:
            return (lons,lats)
        return 

# -------------------------- M A I N -----------------------------
if __name__ == "__main__":
    
    
    # TEMPO based defaults
    # --------------------
    lon_0 = -100.
    bbox = '3112801.8567295657,7356642.8104572548,8731788.567404747,10354334.448304549'
    Nx = 1250 # number of X pixels 
    outDir = '.'
    rPat='_globe_,_tempo_'

#   Parse command line options
#   --------------------------
    parser = OptionParser(usage="Usage: %prog [OPTIONS] geos5_global_image(s)",
                          version='1.0.0' )

    parser.add_option("-b", "--bbox", dest="bbox", default=bbox,
              help="Basemap GEO bounding box in map projection coordinates (default=%s)"\
                          %bbox )

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
              help="Number of NORTH-SOUTH pixels; default is based on Nx and aspect ratio in map projection coordinates" )

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

    # parse bbox
    # ----------
    try:
        bbox = [ float(v) for v in options.bbox.split(',') ]
    except:
        parse.error("invalid bbox")
    if len(bbox) != 4:
        parse.error("bbox must have at leas 4 CSVs")

    # Instantiate basemap transform for GEO
    # -------------------------------------
    m =  MyBasemap(projection='geos',lon_0=lon_0,resolution=None,
                   llcrnrx=bbox[0],
                   llcrnry=bbox[1],
                   urcrnrx=bbox[2],
                   urcrnry=bbox[3],
                   rsphere=rsphere,
                   satellite_height=satellite_height)

    # Define bounding box
    # -------------------
    m.set_bbox(bbox)

    # Cache Coordinates
    # -----------------
    Nx = options.Nx
    if options.Ny == None:
        Ny = int(((bbox[3]-bbox[1])/(bbox[2]-bbox[0]))*Nx)
    else:
        Ny = options.Ny
    Lons, Lats, x, y = m.makegrid(Nx,Ny, returnxy=True, setCache=True)

    # Get coordinates from file (for debugging purposes)
    # --------------------------------------------------
    if False:
        nc = Dataset('/nobackup/TEMPO/LevelG/invariant/tempo.lg1.invariant.nc4')
        Lons = nc.variables['clon'][:,124:1374]
        Lats = nc.variables['clat'][:,124:1374]
        Lons = Lons[-1::-1,-1::-1]
        Lats = Lats[-1::-1,-1::-1]
        m.makegrid_grid = (Lons, Lats, x, y )
        Ny, Nx = Lons.shape

    # Loop over images
    # ----------------
    for imFile in Images:

        # Form output file
        # ----------------
        outFile = options.outDir + '/' + os.path.basename(imFile).replace(p1,p2)

        # Unless force, do not overwrite files
        # ------------------------------------
        if not options.force:
            if os.path.exists(outFile):
                print("[x] Output file exists, skipping %s"%outFile)
                continue

        if options.verbose:
            print("[ ] Reading %s"%os.path.basename(imFile))

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

        # Interpolate from global latlon to GEO sector
        # --------------------------------------------
        rgb = ma.zeros((Ny,Nx,nc))
        for k in range(nc):
            if options.verbose:
                print("   - Remapping Layer ", k)
            rgb[:,:,k] = m.transform_scalar(RGB[-1::-1,:,k],lons,lats,
                                            Nx, Ny,
                                            masked=True)
        # Save the image
        # --------------
        if options.verbose:
            print("   - Writing %s"%outFile)
        imsave(outFile,rgb)
        
