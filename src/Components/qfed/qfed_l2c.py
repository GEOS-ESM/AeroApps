#!/usr/bin/env python

"""
  A Python scrpipt to create QFED Level 2c files.
"""

import warnings
warnings.simplefilter('ignore',DeprecationWarning)

import os
import sys

from numpy     import median, savez
from optparse  import OptionParser   # Command-line args  

from   qfed      import QFED
from   qfed.Date import Date
from   qfed.met  import MET

#---------------------------------------------------------------------

if __name__ == "__main__":

    level2_dir = '/nobackup/2/MODIS/Level2'
    igbp_dir  = '/nobackup/Emissions/Vegetation/GL_IGBP_INPE'
    l2c_file = 'qfed_l2c.txt'
    met_file = '/nobackup/1/MERRA/opendap/chm_Fv_3hourly.ddf'

    products  = 'MOD14,MYD14'
    qc_thresh = 0.    
    pm_top = 50.

#   Parse command line options
#   --------------------------
    parser = OptionParser(usage="Usage: %prog [options] year doy_beg [doy_end]",
                          version='qfed_level2c-1.0.0' )

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
    
    parser.add_option("-o", "--output", dest="l2c_file", default=l2c_file,
                      help="output Level 2C file name (default=%s)"%l2c_file )
    

    parser.add_option("-p", "--products", dest="products", default=products,
                      help="CSV list of MODIS fire products (default=%s)"\
                           %products )

    parser.add_option("-P", "--ptop", dest="ptop", type="float", default=pm_top,
                      help="top pressure in hPa (default=%f)"%pm_top )
    
    parser.add_option("-q", "--qc", dest="qc", type="float", default=qc_thresh,
                      help="q/c threshold (default=%f)"%qc_thresh )
    
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
        print ""
        print "                          QFED Level 2c Processing"
        print "                          ------------------------"
        print ""
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
            fires = QFED(Path,Verb=Verb,qc_thresh=options.qc)
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

#           Instantitel Met fields
#           ---------------------
            met = MET(options.met_file,
                      Verb=Verb, pm_top = options.ptop*100.,
                      Vars=('u','v','t','qv','delp'))

#           Compute plume extent given Met fields
#           -------------------------------------
            fires.getPlume(met)

#           Write out Level 2C file
#           -----------------------
            fires.print_l2c(MxD14+'.'+options.l2c_file)

def save_met():

    savez('met_fields.npz', lon  = met.lon,
                            lat  = met.lat,
                            lev  = met.lev,
                            ak   = met.ak,
                            bk   = met.bk,
                            delp = met.fields['delp'],
                            u    = met.fields['u'],
                            v    = met.fields['v'],
                            T    = met.fields['t'],
                            qv   = met.fields['qv'])

    






