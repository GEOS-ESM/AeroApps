#!/usr/bin/env python3

"""
    Wrapper for stn_sampler.py
    Samples G5NR collections for the data we need to run vlidort simulator

    Patricia Castellanos, May, 2017

"""

import os
import subprocess
import time
import argparse
from datetime        import datetime, timedelta
from dateutil.parser import parse         as isoparser
from MAPL            import Config
from netCDF4         import Dataset
import numpy         as np

if os.path.exists('/discover/nobackup'):
    nccat = '/usr/local/other/SLES11.1/nco/4.4.4/intel-12.1.0.233/bin/ncrcat'
else:
    nccat = '/ford1/share/dasilva/bin/ncrcat'

#------------------------------------ M A I N ------------------------------------

if __name__ == "__main__":

    # Defaults
    DT_hours = 24
    algo     = "linear"

    parser = argparse.ArgumentParser()
    parser.add_argument("iso_t1",
                        help="starting iso time")
    parser.add_argument("iso_t2",
                        help="ending iso time")

    parser.add_argument("prep_config",
                        help="prep config filename")

    parser.add_argument("-a", "--algo", default=algo,
              help="Interpolation algorithm, one of linear, nearest (default=%s)"\
                          %algo)
    parser.add_argument("-D","--DT_hours", default=DT_hours, type=int,
                        help="Timestep in hours for each file (default=%i)"%DT_hours)

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
        rcFiles = rcFiles.replace(' ','').split(',')
    else:
        rcFiles = (rcFiles.replace(' ',''),)

    colNames = cf('COLNAMES')
    if ',' in colNames:
        colNames = colNames.replace(' ','').split(',')
    else:
        colNames = (colNames.replace(' ',''),)

    outdir = cf('OUTDIR')
    instname = cf('INSTNAME')
    stnFile  = cf('STNFILE')
    dt_secs  = cf('DT_SECS')

    date = isoparser(args.iso_t1)
    enddate   = isoparser(args.iso_t2)

    pdt   = timedelta(hours=1)
    while date <= enddate:
        outpath = '{}/Y{}/M{}'.format(outdir,date.year,str(date.month).zfill(2))
        if not os.path.exists(outpath):
            os.makedirs(outpath)

        # run trajectory sampler on model fields
        for rc,colname in zip(rcFiles,colNames):
            nymd = str(date.date()).replace('-','')
            hour = str(date.hour).zfill(2)
            edate = date 
        
            outFile = '{}/{}-g5nr.lb2.{}.{}_{}z.nc4'.format(outpath,instname,colname,nymd,hour)

            Options =     " --output=" + outFile       + \
                          " --isoTime"  +\
                          " --algorithm=" + args.algo  +\
                          " --dt_secs=" + dt_secs

            if args.verbose:
                Options += " --verbose" 

            cmd = 'python -u /discover/nobackup/pcastell/workspace/GAAS/Linux/bin/g5nr_stn_sampler.py {} {} {} {} {}'.format(Options,stnFile,rc,date.isoformat(),edate.isoformat())
            print(cmd)
            if not args.dryrun:
                os.system(cmd)


        date += timedelta(hours=args.DT_hours)

