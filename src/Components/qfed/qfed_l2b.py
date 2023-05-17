#!/usr/bin/env python3

"""
  A Python scrpipt to create QFED Level 2b files.
"""

import warnings
warnings.simplefilter('ignore',DeprecationWarning)

import os
import sys

from numpy     import median, savez
from optparse  import OptionParser   # Command-line args  

import qfed

#---------------------------------------------------------------------

if __name__ == "__main__":

    igbp_dir  = '/nobackup/Emissions/Vegetation/GL_IGBP_INPE'
    l2b_file  = 'qfed_l2b.txt'
    qc_thresh = 0.    
    npz_file = 'NONE'

#   Parse command line options
#   --------------------------
    parser = OptionParser(usage="Usage: %prog [options] MxD14_granule_path",
                          version='qfed_level2b-1.0.0' )

    parser.add_option("-n", "--npz", dest="npz_file", default=npz_file,
                      help="binary numpy npz file name (default=%s)"%npz_file )
    
    parser.add_option("-i", "--igbp", dest="igbp_dir", default=igbp_dir,
                      help="path to IGBP vegetation database (default=%s)"\
                           %igbp_dir )
    
    parser.add_option("-o", "--output", dest="l2b_file", default=l2b_file,
                      help="output Level 2B file name (default=%s)"%l2b_file )
    
    parser.add_option("-q", "--qc", dest="qc", type="float", default=qc_thresh,
                      help="q/c threshold (default=%f)"%qc_thresh )
    
    parser.add_option("-v", "--verbose",
                      action="store_true", dest="verbose")

    (options, mxd14_dir) = parser.parse_args()
    
    if len(mxd14_dir) == 0:
        parser.error("missing granule directory")

    if options.verbose:
        Verb=1
        print("")
        print("                          QFED Level 2b Processing")
        print("")
    else:
        Verb=0

#   Create Level 2B file
#   --------------------
    fires = qfed.QFED(mxd14_dir,Verb=Verb,qc_thresh=options.qc)
    fires.getSimpleVeg(Path=options.igbp_dir)

#   Compute fire size
#   -----------------
    fires.dozier(Verbose=True)  

    m = fires.m

#   Write out Level 2B file
#   -----------------------
    fires.print_l2b(options.l2b_file)

#   Optionally write npz file
#   -------------------------
    if options.npz_file != "NONE":
        savez(options.npz_file, lon   = fires.lon[m],
                                lat   = fires.lat[m],
                                veg   = fires.veg[m],
                                pow   = fires.pow[m],
                                Tf    = fires.Tf[m],
                                farea = fires.farea[m],
                                hflux = fires.hflux[m],
                                pixar = fires.pixar[m])
        print("")
        print("[] wrote "+options.npz_file)


    






