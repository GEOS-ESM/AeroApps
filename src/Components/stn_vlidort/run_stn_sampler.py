#!/usr/bin/env python

"""
    Wrapper for stn_sampler.py
    Samples G5NR collections for the data we need to run vlidort simulator

    Patricia Castellanos, May, 2017

"""

import os
import subprocess
import argparse
from datetime        import datetime, timedelta
from dateutil.parser import parse         as isoparser
from MAPL            import Config

if os.path.exists('/discover/nobackup'):
    nccat = '/usr/local/other/SLES11.1/nco/4.4.4/intel-12.1.0.233/bin/ncrcat'
else:
    nccat = '/ford1/share/dasilva/bin/ncrcat'

#------------------------------------ M A I N ------------------------------------

if __name__ == "__main__":

    # Defaults
    DT_hours = 24
    algo     = "linear"
    nproc    = 4

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

    parser.add_argument("-n", "--nproc",default=nproc,type=int,
                        help="Number of processors (default=%i)."%nproc)    

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

    Date = isoparser(args.iso_t1)
    enddate   = isoparser(args.iso_t2)

    pdt  = timedelta(hours=args.DT_hours/args.nproc)


    while Date < enddate:
        outpath = '{}/Y{}/M{}'.format(outdir,Date.year,str(Date.month).zfill(2))
        if not os.path.exists(outpath):
            os.makedirs(outpath)

        # run trajectory sampler on model fields
        # split across multiple processors by date
        datelist = [Date + p*pdt for p in range(args.nproc)]
        filelist = []
        for rc,colname in zip(rcFiles,colNames):
            processes = set()
            for date in datelist:
                nymd = str(date.date()).replace('-','')
                hour = str(date.hour).zfill(2)
                edate = date + pdt
            
                outFile = '{}/{}-g5nr.lb2.{}.{}_{}z.nc4'.format(outpath,instname,colname,nymd,hour)

                Options =     " --output=" + outFile       + \
                              " --format=NETCDF4_CLASSIC"      + \
                              " --isoTime"  +\
                              " --algorithm=" + args.algo

                if args.verbose:
                    Options += " --verbose" 

                cmd = './g5nr_stn_sampler.py {} {} {} {} {}'.format(Options,stnFile,rc,date.isoformat(),edate.isoformat())
                print cmd
                if not args.dryrun:
                    processes.add(subprocess.Popen(cmd, shell=True))

                filelist.append(outFile)

            #Wait till all the processes are finished
            for p in processes:
                if p.poll() is None:
                    p.wait()

            if not args.dryrun:
                #Concatenate outfiles into one
                cmd = nccat + ' -d time -H -h -c -A ' + ' '.join(filelist) +' -o ' + filelist[0]
                print cmd
                if os.system(cmd):
                    raise ValueError, "nccat failed for {}".format(nymd)


        Date += timedelta(hours=args.DT_hours)

