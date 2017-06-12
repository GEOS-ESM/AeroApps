#!/usr/bin/env python

"""
    Wrapper for stn_sampler.py
    Samples G5NR collections for the data we need to run vlidort simulator
    Special case for invariant collections

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


#------------------------------------ M A I N ------------------------------------
if __name__ == "__main__":

    # Defaults
    DT_hours = 24
    algo     = "linear"
    nproc    = 8

    parser = argparse.ArgumentParser()

    parser.add_argument("prep_config",
                        help="prep config filename")

    parser.add_argument("-a", "--algo", default=algo,
              help="Interpolation algorithm, one of linear, nearest (default=%s)"\
                          %algo)

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

    Date = datetime(2006,01,01,00)
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    # run station sampler on model fields
    for rc,colname in zip(rcFiles,colNames):       
        outFile = '{}/{}-g5nr.lb2.{}.nc4'.format(outdir,instname,colname)

        Options =     " --output=" + outFile       + \
                      " --format=NETCDF4_CLASSIC"      + \
                      " --isoTime"  +     \
                      " --dt_secs=3600"  +\
                      " --algorithm=" + args.algo

        if args.verbose:
            Options += " --verbose" 

        cmd = './g5nr_stn_sampler.py {} {} {} {} {}'.format(Options,stnFile,rc,Date.isoformat(),Date.isoformat())

        if not args.dryrun:
            if os.system(cmd):
                raise ValueError, "g5nr_stn_sampler failed for {}".format(colname)