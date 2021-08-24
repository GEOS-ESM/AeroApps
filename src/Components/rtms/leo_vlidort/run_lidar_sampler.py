#!/usr/bin/env python

"""
    Wrapper for trj_sampler.py
    Samples G5NR collections for the data we need to run vlidort simulator

    Patricia Castellanos, May, 2017

"""

import os
import errno
import subprocess
import argparse
import time
from datetime        import datetime, timedelta
from dateutil.parser import parse         as isoparser
from MAPL            import Config
from netCDF4         import Dataset
import numpy         as np

if os.path.exists('/discover/nobackup'):
    nccat = '/usr/local/other/nco/4.8.1/bin/ncrcat'
else:
    nccat = '/ford1/share/dasilva/bin/ncrcat'

#------------------------------------ M A I N ------------------------------------
def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise

def fix_time(filelist,tbeg):
    for filename in filelist:
        nc = Dataset(filename,'r+')
        time = nc.variables['time']
        time.units = 'seconds since %s'%tbeg.isoformat(' ')
        tyme = nc.variables['isotime'][:]
        tyme = np.array([isoparser(''.join(t)) for t in tyme])
        time[:] = np.array([(t-tbeg).total_seconds() for t in tyme])

        nc.close()


def StartNew(processes,cmds,nextdate,lendate):
   """ Start a new subprocess if there is work to do """

   if nextdate < lendate:
      proc = subprocess.Popen(cmds[nextdate], shell=True)
      print cmds[nextdate]
      nextdate += 1
      processes.append(proc)

   return processes,nextdate

def CheckRunning(processes,cmds,nextdate,lendate,args):
   """ Check any running processes and start new ones if there are spare slots."""

   for p in range(len(processes))[::-1]: # Check the processes in reverse order
      if processes[p].poll() is not None: # If the process hasn't finished will return None
         del processes[p] # Remove from list - this is why we needed reverse order

   while (len(processes) < args.nproc) and (nextdate < lendate): # More to do and some spare slots
      processes, nextdate = StartNew(processes,cmds,nextdate,lendate)

   return processes,nextdate

if __name__ == "__main__":

    # Defaults
    DT_hours = 24
    dt_secs = "1"
    algo    = "linear"
    nproc    = 8
    exp     = 'g5nr'

    parser = argparse.ArgumentParser()
    parser.add_argument("iso_t1",
                        help="starting iso time")
    parser.add_argument("iso_t2",
                        help="ending iso time")

    parser.add_argument("prep_config",
                        help="prep config filename")

    parser.add_argument("-e", "--exp", default=exp,
              help="GEOS experiment name (default=%s)"\
                          %exp)

    parser.add_argument("-a", "--algo", default=algo,
              help="Interpolation algorithm, one of linear, nearest (default=%s)"\
                          %algo)
    parser.add_argument("-D","--DT_hours", default=DT_hours, type=int,
                        help="Timestep in hours for each file (default=%i)"%DT_hours)

    parser.add_argument("-d", "--dt_secs", default=dt_secs,
                        help="Timesetp in seconds for TLE sampling (default=%s)"%dt_secs )    

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
    tleFile  = cf('TLEFILE')

    Date = isoparser(args.iso_t1)
    enddate   = isoparser(args.iso_t2)

    if args.DT_hours == 1:
        #pdt   = timedelta(minutes=int(60/args.nproc))
        pdt   = timedelta(minutes=5)
    else:
        #pdt   = timedelta(hours=int(args.DT_hours/args.nproc))
        pdt   = timedelta(hours=1)

    while Date < enddate:
        if args.DT_hours < 24:
            outpath = '{}/Y{}/M{}/D{}'.format(outdir,Date.year,str(Date.month).zfill(2),str(Date.day).zfill(2))
        else:
            outpath = '{}/Y{}/M{}'.format(outdir,Date.year,str(Date.month).zfill(2))
        mkdir_p(outpath)

        # run trajectory sampler on model fields
        # split across multiple processors by date
        if args.DT_hours == 1:
                datelist = [Date + p*pdt for p in range(60/5)]
        else:        
        	datelist = [Date + p*pdt for p in range(args.DT_hours)]
        lendate  = len(datelist)

        for rc,colname in zip(rcFiles,colNames):
            processes = []
            cmds      = []
            filelist  = []
            for date in datelist:
                nymd   = str(date.date()).replace('-','')
                hour   = str(date.hour).zfill(2)
                minute = str(date.minute).zfill(2)

                edate = date + pdt - timedelta(seconds=int(args.dt_secs))

                outFile = '{}/{}-{}.lb2.{}.{}_{}{}z.nc4'.format(outpath,instname,args.exp,colname,nymd,hour,minute)

                Options =     " --rcFile=" + rc      + \
                              " --output=" + outFile       + \
                              " --format=NETCDF4_CLASSIC"      + \
                              " --dt_secs=" + args.dt_secs      + \
                              " --isoTime"  +\
                              " --trajectory=tle"  +\
                              " --algorithm=" + args.algo

                if args.verbose:
                    Options += " --verbose" 

                cmd = './lidar_sampler.py {} {} {} {}'.format(Options,tleFile,date.isoformat(),edate.isoformat())
                cmds.append(cmd)
                filelist.append(outFile)

            # Manage processes
            # This will start the max processes running    
            processes, nextdate = CheckRunning(processes,cmds,0,lendate,args)
            while len(processes)>0: # Some things still going on
                time.sleep(10)      # Wait
                # add more processes as other ones finish
                processes, nextdate = CheckRunning(processes,cmds,nextdate,lendate,args)

            if (not args.dryrun) & (args.nproc > 1):
                # make time units all the same
                fix_time(filelist,Date)
                #Concatenate outfiles into one
                cmd = nccat + ' -h -A ' + ' '.join(filelist) +' -o ' + filelist[0]
                print cmd
                if os.system(cmd):
                    raise ValueError, "nccat failed for {}".format(nymd)

                for filename in filelist[1:]:
                    os.remove(filename)

        Date += timedelta(hours=args.DT_hours)

