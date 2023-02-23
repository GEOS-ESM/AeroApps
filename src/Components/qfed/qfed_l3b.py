#!/usr/bin/env python3
# -W ignore::DeprecationWarning

"""
  A Python script to create QFED Level 3b files.
"""

import warnings
warnings.simplefilter('ignore',DeprecationWarning)
warnings.simplefilter('always',UserWarning)

import os
from datetime       import date, timedelta

from optparse       import OptionParser
from glob           import glob

from gfio           import GFIO
from qfed.emissions import Emissions

Sat = { 'MOD14': 'MODIS_TERRA', 'MYD14': 'MODIS_AQUA' }

#---------------------------------------------------------------------

if __name__ == "__main__":

    expid       = 'qfed2'
    level3a_dir = '/nobackup/2/MODIS/Level3'
    level3b_dir = '/nobackup/2/MODIS/Level3'
    products    = 'MOD14,MYD14'
    fill_days   = 1

    
#   Parse command line options
#   --------------------------
    parser = OptionParser(usage="Usage: %prog [options] year doy_beg [doy_end]",
                          version='qfed_level3b-1.0.b2' )

    parser.add_option("-i", "--input", dest="level3a_dir",
                      default=level3a_dir,
                      help="directory for input FRP files (default=%s)"\
                           %level3a_dir )

    parser.add_option("-o", "--output", dest="level3b_dir",
                      default=level3b_dir,
                      help="directory for output Emission files (default=%s)"\
                           %level3b_dir )

    parser.add_option("-p", "--products", dest="products", default=products,
                      help="CSV list of MODIS fire products (default=%s)"\
                           %products )
    
    parser.add_option("-x", "--expid", dest="expid", default=expid,
                      help="experiment id (default=%s)"\
                           %expid )

    parser.add_option("-u", "--uncompressed",
                      action="store_true", dest="uncompressed",
                      help="do not use n4zip to compress gridded output file (default=False)")

    parser.add_option("-n", "--ndays", dest="ndays", type="int", default=fill_days,
                      help="Number of days to fill in (default=%d)"\
                           %fill_days )

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
        print("                          QFED Level 3B Processing")
        print("                          ------------------------")
        print("")
    else:
        Verb=0

    Products = options.products.split(',')

#   Grid FRP and observed area
#   --------------------------
    for doy in range(int(doy_beg),int(doy_end)+1):

        d = date(int(year),1,1) + timedelta(days=(doy - 1))

#       Read input FRP and Area for each satellite
#       ------------------------------------------
        FRP, Land, Water, Cloud, F = ( {}, {}, {}, {}, {} )
        for MxD14 in Products:

            sat = Sat[MxD14]

#           Open input file
#           ---------------
            l3a_dir  = os.path.join(options.level3a_dir, MxD14, 'Y%04d'%d.year, 'M%02d'%d.month)
            l3a_file = '%s_%s.frp.???.%04d%02d%02d.nc4'%(options.expid, MxD14, d.year, d.month, d.day)
            
            pat = os.path.join(l3a_dir, l3a_file)

            try:
                ifn = glob(pat)[0]
                f = GFIO(ifn)
            except:
                print("[x] cannot find/read input FRP file for %s, ignoring it"%d)
                continue

            if Verb:
                print("[] Reading ", ifn) 

            Land[sat]  = f.read('land')
            Water[sat] = f.read('water')
            Cloud[sat] = f.read('cloud')
            
            FRP[sat] = [ f.read('frp_tf'),
                         f.read('frp_xf'),
                         f.read('frp_sv'),
                         f.read('frp_gl') ]

            F[sat]   = [ f.read('fb_tf'),
                         f.read('fb_xf'),
                         f.read('fb_sv'),
                         f.read('fb_gl') ]

            col = ifn.split('/')[-1].split('.')[2] # collection


#       FRP density forecast files
#       --------------------------
        d_ = d + timedelta(days=1)
        forecast_files = {}
        for MxD14 in Products:
            sat = Sat[MxD14]
 
            if sat in list(FRP.keys()):
                l3a_dir  = os.path.join(options.level3a_dir, MxD14, 'Y%04d'%d_.year, 'M%02d'%d_.month)
                l3a_file = '%s_%s.frp.%s.%04d%02d%02d.nc4'%(options.expid, MxD14, col, d_.year, d_.month, d_.day)
            
                forecast_files[sat] = os.path.join(l3a_dir, l3a_file)


#       Create the top level directory for output files
#       -----------------------------------------------
        dir = os.path.join(options.level3b_dir, 'QFED')
        rc = os.system("/bin/mkdir -p %s"%dir)
        if rc:
            raise IOError('cannot create output directory')

#       Write output file
#       -----------------
        E = Emissions(d, FRP, F, Land, Water, Cloud, Verb=Verb)
        E.calculate()
        E.write(dir=dir, forecast=forecast_files, expid=options.expid, col=col, ndays=options.ndays, 
                uncompressed=options.uncompressed)
        
