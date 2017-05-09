#!/usr/bin/env python

"""
    Wrapper for trj_sampler.py
    Samples G5NR collections for the data we need to run vlidort simulator

    Patricia Castellanos, May, 2017

"""

import os
import argparse
from datetime        import datetime, timedelta
from dateutil.parser import parse         as isoparser
from MAPL            import Config

#------------------------------------ M A I N ------------------------------------

if __name__ == "__main__":

    # Defaults
    DT_hours = 24
    dt_secs = "60"

    parser = argparse.ArgumentParser()
    parser.add_argument("iso_t1",
                        help="starting iso time")
    parser.add_argument("iso_t2",
                        help="ending iso time")

    parser.add_argument("prep_config",
                        help="prep config filename")


    parser.add_argument("-D","--DT_hours", default=DT_hours, type=int,
                        help="Timestep in hours for each file (default=%i)"%DT_hours)

    parser.add_argument("-d", "--dt_secs", default=dt_secs,
                        help="Timesetp in seconds for TLE sampling (default=%s)"%dt_secs )    

    parser.add_argument("-v", "--verbose",action="store_true",
                        help="Verbose mode (default=False).")

    parser.add_argument("-r", "--dryrun",action="store_true",
                        help="do a dry run (default=False).")    



    args = parser.parse_args()

    # Parse prep config
    # -----------------
    cf = Config(args.prep_config,delim=' = ')

    rcFiles = cf('RCFILES')
    if ',' in rcFiles:
        rcFiles = rcFiles.split(',')
    else:
        rcFiles = (rcFiles,)


    outdir = cf('OUTDIR')
    instname = cf('INSTNAME')
    tleFile  = cf('TLEFILE')

    date = isoparser(args.iso_t1)
    enddate   = isoparser(args.iso_t2)

    while date < enddate:
        nymd = str(date.date()).replace('-','')
        hour = str(date.hour).zfill(2)
        outpath = '{}/Y{}/M{}'.format(outdir,date.year,str(date.month).zfill(2))
        if not os.path.exists(outpath):
            os.makedirs(outpath)
        for rc in rcFiles:
            colname = '_'.join(rc.split('_')[0:2])

            outFile = '{}/{}-g5nr.lb2.{}.{}_{}z.nc4'.format(outpath,instname,colname,nymd,hour)

            Options =     " --rcFile=" + rc      + \
                          " --output=" + outFile       + \
                          " --format=NETCDF4_CLASSIC"      + \
                          " --dt_secs=" + args.dt_secs      + \
                          " --isoTime"  +\
                          " --trajectory=tle"

            if args.verbose:
                Options += " --verbose" 


            cmd = 'trj_sampler.py {} {} {} {}'.format(Options,tleFile,date.isoformat(),(date+timedelta(hours=args.DT_hours-1)).isoformat())
            print cmd
            if not args.dryrun:
                if os.system(cmd):
                    raise ValueError, "trj_sampler.py failed for %s on %s"%(rc, date)            


        date += timedelta(hours=args.DT_hours)

