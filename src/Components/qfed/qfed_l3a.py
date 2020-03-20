#!/usr/bin/env python
# -W ignore::DeprecationWarning

"""
  A Python script to create QFED Level 3a files.
"""

import warnings
warnings.simplefilter('ignore',DeprecationWarning)
warnings.simplefilter('always',UserWarning)

import os

from datetime      import date, timedelta

from optparse      import OptionParser   # Command-line args  

from qfed.mxd14_l3 import MxD14_L3


#---------------------------------------------------------------------

if __name__ == "__main__":

    expid      = 'qfed2'
    level1_dir = '/nobackup/2/MODIS/Level1'
    level2_dir = '/nobackup/2/MODIS/Level2'
    level3_dir = '/nobackup/2/MODIS/Level3'
    products   = 'MOD14,MYD14'
    igbp_dir   = '/nobackup/Emissions/Vegetation/GL_IGBP_INPE'
    res        = 'e'
    
#   Parse command line options
#   --------------------------
    parser = OptionParser(usage="Usage: %prog [options] year doy_beg [doy_end]",
                          version='qfed_level3a-1.0.b2' )

    parser.add_option("-f", "--mxd14", dest="level2_dir", default=level2_dir,
                      help="directory for MOD14/MYD14 fire files (default=%s)"\
                           %level2_dir )

    parser.add_option("-g", "--mxd03", dest="level1_dir", default=level1_dir,
                      help="directory for MOD03/MYD03 geolocaltion files (default=%s)"\
                           %level1_dir )

    parser.add_option("-i", "--igbp", dest="igbp_dir", default=igbp_dir,
                      help="directory for IGBP vegetation database (default=%s)"\
                           %igbp_dir )
    
    parser.add_option("-o", "--output", dest="level3_dir", default=level3_dir,
                      help="directory for output files (default=%s)"\
                           %level3_dir )

    parser.add_option("-p", "--products", dest="products", default=products,
                      help="CSV list of MODIS fire products (default=%s)"\
                           %products )
    
    parser.add_option("-r", "--resolution", dest="res", default=res,
                      help="horizontal resolution: a for 4x5, b for 2x2.5, etc. (default=%s)"\
                           %res )
    
    parser.add_option("-x", "--expid", dest="expid", default=expid,
                      help="experiment id (default=%s)"\
                           %expid )

    parser.add_option("-q", "--disable-qc", 
                      action="store_true", dest="disable_qc",
                      help="disable quality control procedures (default=%s)"\
                           % False)

    parser.add_option("-b", "--bootstrap", 
                      action="store_true", dest="bootstrap",
                      help="initialize FRP forecast (default=False)")

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
        print ""
        print "                          QFED Level 3A Processing"
        print "                          ------------------------"
        print ""
    else:
        Verb=0

    Products = options.products.split(',')

#   Grid FRP and observed area
#   --------------------------
    for doy in range(int(doy_beg), int(doy_end)+1):
        d = date(int(year), 1, 1) + timedelta(days=(doy - 1))
        d_= d + timedelta(days=1)

#       Loop over products
#       ------------------
        for MxD14 in Products:

            Path = os.path.join(options.level2_dir, MxD14, '%04d'%d.year, '%03d'%doy)

#           Do the gridding for this product
#           --------------------------------
            fires = MxD14_L3(Path,
                             options.level1_dir,
                             options.igbp_dir,
                             res=options.res,
                             Verb=Verb)

#           Create directory for output file
#           --------------------------------
            dir_a = os.path.join(options.level3_dir, MxD14, 'Y%04d'%d.year, 'M%02d'%d.month)
            dir_f = os.path.join(options.level3_dir, MxD14, 'Y%04d'%d_.year, 'M%02d'%d_.month)

            dir = {'ana': dir_a, 'bkg': dir_f}
            for k in dir.keys():
                rc  = os.system("/bin/mkdir -p %s"%dir[k])
                if rc:
                    raise IOError, "Cannot create output directory '%s'" % dir[k]


#           Quality Control
#           ---------------
            qc = not options.disable_qc

#           Write output file
#           -----------------
            fires.write(dir=dir, expid=options.expid+'_'+MxD14, qc=qc, bootstrap=options.bootstrap)

#           Compress it unless the user disabled it
#           ---------------------------------------
            if not options.uncompressed and fires.doy != None:
                rc = os.system("n4zip %s"%fires.filename)
                if rc:
                    warnings.warn('cannot compress output file <%s>'%fires.filename)

