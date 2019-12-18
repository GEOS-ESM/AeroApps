#!/usr/bin/env python

"""
    Wrapper for mcd12c.py
    Samples MCD12C1 collection on a lidar trajectory.
    Lidar_sampler.py has already been called

    Patricia Castellanos, May, 2017

"""

import os
import argparse
from datetime        import datetime, timedelta
from dateutil.parser import parse         as isoparser
from MAPL            import Config
from mcd12c          import MCD12C

#------------------------------------ M A I N ------------------------------------

if __name__ == "__main__":
    # Defaults
    DT_hours = 24

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



    args = parser.parse_args()

    # Parse prep config
    # -----------------
    cf = Config(args.prep_config,delim=' = ')
    outdir    = cf('OUTDIR')
    indir     = cf('INDIR')
    datadir   = cf('DATADIR')
    instname  = cf('INSTNAME')

    date = isoparser(args.iso_t1)
    enddate   = isoparser(args.iso_t2)

    while date < enddate:
        nymd = str(date.date()).replace('-','')
        hour = str(date.hour).zfill(2)
        outpath = '{}/Y{}/M{}'.format(outdir,date.year,str(date.month).zfill(2))
        inpath  = '{}/Y{}/M{}'.format(indir,date.year,str(date.month).zfill(2))
        if not os.path.exists(outpath):
            os.makedirs(outpath)


        outFile = '{}/{}-g5nr.lb2.land_cover.{}_{}z.nc4'.format(outpath,instname,nymd,hour)
        inFile  = '{}/{}-g5nr.lb2.aer_Nv.{}_{}z.nc4'.format(inpath,instname,nymd,hour)

        # Read in trajectory lat/lon
        # Interpolate BRDF 

        brdf = MCD12C(datadir)
        brdf.sample(inFile,outFile,Verbose=args.verbose)

        date += timedelta(hours=args.DT_hours)


