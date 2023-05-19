#!/usr/bin/env python3

"""
  A Python scrpipt to create QFED Level 3C files.
"""

import warnings
warnings.simplefilter('ignore',DeprecationWarning)

import os
import sys

from numpy     import median, savez
from optparse  import OptionParser   # Command-line args  

from   qfed      import PLUME_L2, PLUME_L3
from   qfed.met  import MET
from   qfed.Date import Date

#---------------------------------------------------------------------

if __name__ == "__main__":

    expid = 'qfed2'
    level2_dir = '/nobackup/2/MODIS/Level2'
    level3_dir = '/nobackup/2/MODIS/Level3'
    igbp_dir   = '/nobackup/Emissions/Vegetation/GL_IGBP_INPE'
    met_file   = '/nobackup/1/MERRA/opendap/chm_Fv_3hourly.ddf'

    products   = 'MOD14,MYD14'
    qc_thresh  = 0.    
    pm_top = 50.
    res = 'e'

#   Parse command line options
#   --------------------------
    parser = OptionParser(usage="Usage: %prog [options] year doy_beg [doy_end]",
                          version='qfed_level3c-1.0.0' )

    parser.add_option("-D", "--debug",
                      help="debug mode, process only 99 fires",
                      action="store_true", dest="debug")

    parser.add_option("-f", "--mxd14", dest="level2_dir", default=level2_dir,
                      help="directory for MOD14/MYD14 fire files (default=%s)"\
                           %level2_dir )

    parser.add_option("-i", "--igbp", dest="igbp_dir", default=igbp_dir,
                      help="path to IGBP vegetation database (default=%s)"\
                           %igbp_dir )
    
    parser.add_option("-m", "--met", dest="met_file", default=met_file,
                      help="input Met file name (default=%s)"%met_file )
    
    parser.add_option("-o", "--output", dest="level3_dir", default=level3_dir,
                      help="directory for output files (default=%s)"\
                           %level3_dir )

    parser.add_option("-p", "--products", dest="products", default=products,
                      help="CSV list of MODIS fire products (default=%s)"\
                           %products )

    parser.add_option("-P", "--ptop", dest="ptop", type="float", default=pm_top,
                      help="top pressure in hPa (default=%f)"%pm_top )
    
    parser.add_option("-r", "--resolution", dest="res", default=res,
                      help="horizontal resolution: a for 4x5, b for 2x2.5, etc. (default=%s)"\
                           %res )
    
    parser.add_option("-q", "--qc", dest="qc", type="float", default=qc_thresh,
                      help="q/c threshold (default=%f)"%qc_thresh )
    
    parser.add_option("-x", "--expid", dest="expid", default=expid,
                      help="experiment id (default=%s)"\
                           %expid )

    parser.add_option("-u", "--uncompressed",
                      action="store_true", dest="uncompressed",
                      help="do not use n4zip to compress gridded output file (default=False)")

    parser.add_option("-v", "--verbose",
                      action="store_true", dest="verbose")

    (options, args) = parser.parse_args()
    
    if len(args) == 2:
        year, doy_beg = args
        doy_end = doy_beg
    elif len(args) == 3:
        year, doy_beg, doy_end = args
    else:
        parser.error("must have 2 or 3 arguments: year and day-of-year range")
        
    if options.verbose:
        Verb=1
        print("")
        print("                          QFED Level 3C Processing")
        print("                          ------------------------")
        print("")
    else:
        Verb=0

    Products = options.products.split(',')

#   Loop over days
#   --------------
    for doy in range(int(doy_beg),int(doy_end)+1):

        date = Date((int(year),1,1)) + doy - 1
        doy = "%03d"%doy

#       Loop over products
#       ------------------
        for MxD14 in Products:

            Path = [ level2_dir+'/'+MxD14+'/'+year+'/'+doy,]

#           Read in fire properties, attach biome information
#           -------------------------------------------------
            fires = PLUME_L2(Path,Verb=Verb,qc_thresh=options.qc)
            fires.getSimpleVeg(Path=options.igbp_dir)

#           Compute fire size
#           -----------------
            fires.dozier(Verbose=options.verbose)  

#           Replace fire size with median size for those
#            fires that did not converge
#           --------------------------------------------
            m = fires.m
            i = (m == False)
            fires.p[i] = median(fires.p[m])

#           Short cut for debugging
#           -----------------------
            if options.debug:
                fires.lon = fires.lon[0:100]
                fires.lat = fires.lat[0:100]
                fires.veg = fires.veg[0:100]
                fires.pow = fires.pow[0:100]

#           Instantitel Met fields
#           ---------------------
            met = MET(options.met_file,
                      Verb=Verb, pm_top = options.ptop*100.,
                      Vars=('u','v','t','qv','delp'))

#           Grid fire properties
#           --------------------
            Plume = PLUME_L3(fires,res=options.res)

#           Compute GRIDDED plume extent given Met fields
#           ---------------------------------------------
            Plume.getPlume(met)

#           Create directory for output file
#           --------------------------------
            dir = options.level3_dir+'/'+MxD14+'/Y'+year+'/M%02d'%date.month
            rc = os.system("/bin/mkdir -p %s"%dir)
            if rc:
                raise IOError('cannot create output directory')
            
#           Write output file
#           -----------------
            Plume.write(dir=dir,expid=options.expid+'_'+MxD14)

#           Compress it unless the user disabled it
#           ---------------------------------------
            if not options.uncompressed:
                rc = os.system("n4zip %s"%Plume.filename)
                if rc:
                    warnings.warn('cannot compress output file <%s>'%Plume.filename)





