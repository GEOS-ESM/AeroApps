#!/usr/bin/env python

"""
    Utility to create GEO time dependent geolocation file given an
    invariant file.
"""

import os

from numpy import zeros, ones, arange, array

from optparse        import OptionParser
from datetime        import datetime, timedelta
from dateutil.parser import parse         as isoparser

from netCDF4 import Dataset

def gteInvariant(inFile):
    """
    Read invariant geolocation file and return 
       clon
       clat
       sat_lon
       sat_lat
       sat_alt
    """
    nc = Dataset(inFile)
    
   


# -------------------------- M A I N -----------------------------
if __name__ == "__main__":
    
    inFile = '/nobackup/TEMPO/LG/invariant/tempo.lg1.invariant.nc4'
    outTmpl = 'tempo.lg2.%y4%m2%d2_%h2z.nc4'
    invariant = 'tempo.lg1.invariant.nc4'
    
#   Parse command line options
#   --------------------------
    parser = OptionParser(usage="Usage: %prog [OPTIONS] isoTime",
                          version='1.0.0' )

    parser.add_option("-i", "--invariant", dest="inFile", default=inFile,
              help="Input Invariant LG1 NetCDF file (default=%s)"\
                          %inFile )

    parser.add_option("-v", "--verbose",
                      action="store_true", dest="verbose",
                      help="Verbose mode.")

    # Parse arguments
    # ---------------
    (options, args) = parser.parse_args()
    if len(args) == 1:
        isoTime = args[0]
        tyme = isoparser(isoTime)
    else:
        parser.error("must have 1 arguments: isoTime")

    # Form output file name
    # ---------------------
    outFile = outTmpl.replace('%y4',str(tyme.year)) \
                     .replace('%m2','%02d'%tyme.month) \
                     .replace('%d2','%02d'%tyme.day) \
                     .replace('%h2','%02d'%tyme.hour)

    # Get invariant coordinates
    # -------------------------
    clon, clat, sat_lon, sat_lat, sat_alt = getInvariant(inFile)

    

