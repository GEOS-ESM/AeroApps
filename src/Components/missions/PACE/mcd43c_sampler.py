#!/usr/bin/env python3

"""
    Wrapper for mcd43c.py
    Samples MCD43C1 collection on a lidar trajectory.
    Lidar_sampler.py has already been called

    Patricia Castellanos, May, 2017

"""

import os
import argparse
from datetime        import datetime, timedelta
from dateutil.parser import parse         as isoparser
from MAPL            import Config
from mcd43c          import MCD43C

#------------------------------------ M A I N ------------------------------------

if __name__ == "__main__":
    # Defaults
    DT_mins = 5

    parser = argparse.ArgumentParser()
    parser.add_argument("iso_t1",
                        help="starting iso time")
    parser.add_argument("iso_t2",
                        help="ending iso time")

    parser.add_argument("prep_config",
                        help="prep config filename")

    parser.add_argument("-D","--DT_mins", default=DT_mins, type=int,
                        help="Timestep in minutes for each file (default=%i)"%DT_mins)

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

    while date <= enddate:
        nymd = str(date.date()).replace('-','')
        hms  = date.strftime('%H%M%S')
        outpath = '{}/Y{}/M{}/D{}'.format(outdir,date.year,str(date.month).zfill(2),str(date.day).zfill(2))
        inpath  = '{}/Y{}/M{}/D{}'.format(indir,date.year,str(date.month).zfill(2),str(date.day).zfill(2))
        if not os.path.exists(outpath):
            os.makedirs(outpath)


        outFile = '{}/{}-g5nr.lb.brdf.{}_{}.nc4'.format(outpath,instname,nymd,hms)
        inFile  = '{}/{}-g5nr.lb.asm_Nx.{}_{}.nc4'.format(inpath,instname,nymd,hms)

        # Read in trajectory lat/lon
        # Interpolate BRDF 

        brdf = MCD43C(datadir)
        brdf.sample(inFile,outFile,Verbose=args.verbose)

        date += timedelta(minutes=args.DT_mins)


