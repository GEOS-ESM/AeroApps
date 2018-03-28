#!/usr/bin/env python

"""
    Wrapper for polar_vlidort.py
    Calculated aerosol size distribution

    Patricia Castellanos, March, 2018

"""

import os
import argparse
from   datetime        import datetime, timedelta
from   dateutil.parser import parse         as isoparser
from   MAPL            import Config
from   polar_vlidort   import POLAR_VLIDORT, get_chd
import numpy  as np

#------------------------------------ M A I N ------------------------------------

if __name__ == "__main__":

    # Defaults
    DT_hours = 24

#   Parse command line options
#   --------------------------

    parser = argparse.ArgumentParser()
    parser.add_argument("iso_t1",
                        help="starting iso time")
    parser.add_argument("iso_t2",
                        help="ending iso time")

    parser.add_argument("prep_config",
                        help="prep config filename")

    parser.add_argument("-D","--DT_hours", default=DT_hours, type=int,
                        help="Timestep in hours for each file (default=%i)"%DT_hours)

    parser.add_argument("-v", "--verbose",action="store_true",
                        help="Verbose mode (default=False).")

    parser.add_argument("-r", "--dryrun",action="store_true",
                        help="do a dry run (default=False).")    

    args = parser.parse_args()

    # Parse prep config
    # -----------------
    cf             = Config(args.prep_config,delim=' = ')
    inTemplate     = cf('inDir')     + '/' + cf('inFile')         
    outTemplate    = cf('outDir')    + '/' + cf('outFile')
    channel        = 470
    rcFile         = 'Aod_EOS.rc'

    VZA = None
    HGT = None



    # Change wavelength number in Aod_EOS.rc
    # ---------------------------------------
    channel = '470'
    rcFile = 'Aod_EOS.rc'

    # Loop through dates, running VLIDORT
    # ------------------------------------
    date      = isoparser(args.iso_t1)
    enddate   = isoparser(args.iso_t2)
    Dt        = timedelta(hours=args.DT_hours)

    extOnly  = False
    distOnly = True

    while date < enddate:
        nymd = str(date.date()).replace('-','')
        year = str(date.year)
        month = str(date.month).zfill(2)
        hour = str(date.hour).zfill(2)    

        inFile      = inTemplate.replace('%year',year).replace('%month',month).replace('%nymd',nymd).replace('%hour',hour)
        outFile     = None
        outFileDist = outTemplate.replace('%year',year).replace('%month',month).replace('%nymd',nymd).replace('%hour',hour)

        brdfFile    = None
        ndviFile    = None
        lcFile      = None
        lerFile     = None

        # Initialize VLIDORT class getting aerosol optical properties
        # -----------------------------------------------------------
        print '++++Running VLIDORT with the following arguments+++'
        print '>>>inFile:    ',inFile
        print '>>>outFile:   ',outFile
        print '>>>rcFile:    ',rcFile
        print '>>>albedoType:','None'
        print '>>>channel:   ',channel
        print '>>>VZA:       ',VZA
        print '>>>HGT:       ',HGT
        print '>>>brdfFile:  ',brdfFile
        print '>>>ndviFile:  ',ndviFile
        print '>>>lcFile:    ',lcFile
        print '>>>lerFile    ',lerFile
        print '>>>verbose:   ',args.verbose
        print '++++End of arguments+++'
        if not args.dryrun:
            vlidort = POLAR_VLIDORT(inFile,outFile,rcFile,
                                    None,
                                    channel,
                                    VZA,
                                    HGT,
                                    brdfFile=brdfFile,
                                    ndviFile=ndviFile,
                                    lcFile=lcFile,
                                    lerFile=lerFile,
                                    verbose=args.verbose,
                                    extOnly=extOnly,
                                    distOnly=distOnly,
                                    outFileDist=outFileDist)

            # Run size distribution calculator
            vliort.sizeDistribution()
            vlidort.writeNCdist()

        date += Dt
