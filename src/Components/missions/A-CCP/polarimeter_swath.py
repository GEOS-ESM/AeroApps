#!/usr/bin/env python

"""
    Puts all needed lidar signal simulator variables in one file

    Patricia Castellanos, Dec, 2019

"""

import os
import argparse
from   datetime        import datetime, timedelta
from   dateutil.parser import parse         as isoparser
from   MAPL            import Config
import numpy  as np











#------------------------------------ M A I N ------------------------------------

if __name__ == "__main__":

    # Defaults
    DT_hours = 1

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
    instname       = cf('instname')
    INSTNAME       = instname.upper()
    HGT            = float(cf('HGT'))
    polarimeter_name       = cf('polarimeter_name')


    # Loop through dates, calculating geometry
    # ------------------------------------
    date      = isoparser(args.iso_t1)
    enddate   = isoparser(args.iso_t2)
    Dt        = timedelta(hours=args.DT_hours)

    while date < enddate:
        nymd  = str(date.date()).replace('-','')
        year  = str(date.year)
        month = str(date.month).zfill(2)
        day   = str(date.day).zfill(2)
        hour  = str(date.hour).zfill(2)    

        inFile     = inTemplate.replace('%col','aer_Nv').replace('%year',year).replace('%month',month).replace('%day',day).replace('%nymd',nymd).replace('%hour',hour).replace('%instname',instname).replace('%INSTNAME',INSTNAME)
        outFile    = inTemplate.replace('%col',polarimeter_name).replace('%year',year).replace('%month',month).replace('%day',day).replace('%nymd',nymd).replace('%hour',hour).replace('%instname',instname).replace('%INSTNAME',INSTNAME)


        # Initialize MERGE class and write new outfile
        # -----------------------------------------------------------
        print '++++Getting polarimeter swath with the following arguments+++'
        print '>>>inFile:    ',inFile
        print '>>>outFile:   ',outFile
        print '>>>HGT:       ',HGT
        print '>>>verbose:   ',args.verbose
        print '++++End of arguments+++'
        if not args.dryrun:
            swath = SWATH(inFile,