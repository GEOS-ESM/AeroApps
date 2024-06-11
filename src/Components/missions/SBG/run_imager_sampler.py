#!/usr/bin/env python3

"""
    Wrapper for imager_sampler.py
    Samples G5NR collections for the data we need to run vlidort simulator

    Patricia Castellanos, Jun 2024
    Patricia Castellanos, May 2017

"""

import os
import errno
import subprocess
import argparse
import time
from datetime        import datetime, timedelta
from dateutil.parser import parse         as isoparser
from MAPL.config     import Config
from netCDF4         import Dataset
import numpy         as np

if os.path.exists('/discover/nobackup'):
    nccat = '/usr/local/other/nco/4.8.1/bin/ncrcat'
else:
    nccat = '/ford1/share/dasilva/bin/ncrcat'

def StartNew(processes,cmds,nextproc,totalproc):
   """ Start a new subprocess if there is work to do """

   if nextproc < totalproc:
      proc = subprocess.Popen(cmds[nextproc], shell=True)
      print(cmds[nextproc])
      nextproc += 1
      processes.append(proc)

   return processes,nextproc

def CheckRunning(processes,cmds,nextproc,totalproc,args):
   """ Check any running processes and start new ones if there are spare slots."""

   for p in range(len(processes))[::-1]: # Check the processes in reverse order
      if processes[p].poll() is not None: # If the process hasn't finished will return None
         del processes[p] # Remove from list - this is why we needed reverse order

   while (len(processes) < args.nproc) and (nextproc < totalproc): # More to do and some spare slots
      processes, nextdate = StartNew(processes,cmds,nextproc,totalproc)

   return processes,nextdate
#------------------------------------ M A I N ------------------------------------
if __name__ == "__main__":

    # Defaults
    DT_mins = 5
    algo    = "linear"
    nproc    = 120
    exp     = 'g5nr'

    parser = argparse.ArgumentParser()
    parser.add_argument("iso_t1",
                        help="starting iso time")
    parser.add_argument("iso_t2",
                        help="ending iso time")

    parser.add_argument("track_pcf",
                        help="prep config file with track input file names")

    parser.add_argument("sampler_pcf",
                        help="prep config file with collections to sample")

    parser.add_argument("-e", "--exp", default=exp,
              help="GEOS experiment name (default=%s)"\
                          %exp)

    parser.add_argument("-a", "--algo", default=algo,
              help="Interpolation algorithm, one of linear, nearest (default=%s)"\
                          %algo)
    parser.add_argument("-D","--DT_mins", default=DT_mins, type=int,
                        help="Timestep in minutes for each granule file (default=%i)"%DT_mins)

    parser.add_argument("-v", "--verbose",action="store_true",
                        help="Verbose mode (default=False).")

    parser.add_argument("-r", "--dryrun",action="store_true",
                        help="do a dry run (default=False).")    

    parser.add_argument("-n", "--nproc",default=nproc,type=int,
                        help="Number of processors (default=%i)."%nproc)    

    args = parser.parse_args()

    # Parse prep config
    # -----------------
    cf = Config(args.sampler_config,delim=' = ')

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

    orbitname = cf('ORBITNAME')
    instname  = cf('INSTNAME')

    cf             = Config(args.track_pcf,delim=' = ')
    inTemplate     = cf('inDir')     + '/' + cf('inFile')
    inTemplate     = inTemplate.replace('%orbitname',orbitname).replace('%ORBITNAME',ORBITNAME)

    cmds      = []
    filelist  = []

    Date = isoparser(args.iso_t1)
    enddate   = isoparser(args.iso_t2)
    dt   = timedelta(minutes=args.DT_mins)

    # make a list of proccesses to run
    while Date < enddate:
        nymd   = str(date.date()).replace('-','')
        year  = str(date.year)
        month = str(date.month).zfill(2)
        day   = str(date.day).zfill(2)
        hour   = str(date.hour).zfill(2)
        minute = str(date.minute).zfill(2)

        swathFile = inTemplate.replace('%col',instname).replace('%year',year).replace('%month',month).replace('%day',day).replace('%nymd',nymd).replace('%hour',hour).replace('%minute',minute)
        for rc,colname in zip(rcFiles,colNames):
            nymd   = str(date.date()).replace('-','')
            year  = str(date.year)
            month = str(date.month).zfill(2)
            day   = str(date.day).zfill(2)
            hour   = str(date.hour).zfill(2)
            minute = str(date.minute).zfill(2)

            outFile = inTemplate.replace('%col',colname).replace('%year',year).replace('%month',month).replace('%day',day).replace('%nymd',nymd).replace('%hour',hour).replace('%minute',minute)

            Options =     " --outFile=" + outFile       + \
                          " --format=NETCDF4_CLASSIC"      + \
                          " --isoTime"  +\
                          " --algorithm=" + args.algo

            if args.verbose:
                Options += " --verbose"

            cmd = './imager_sampler.py {} {} {}'.format(Options,swathFile,rc)
            cmds.append(cmd)
            filelist.append(outFile)

        Date += dt

        totalproc = len(cmds)

        # run trajectory sampler on model fields
        # split across multiple processors 
        # This next bit of code manages the processes
        # This will start the max processes running  
        processes = [] 
        processes, nextproc = CheckRunning(processes,cmds,0,totalproc,args)
        while len(processes)>0: # Some things still going on
            time.sleep(10)      # Wait
            # add more processes as other ones finish
            processes, nextproc = CheckRunning(processes,cmds,nextproc,totalproc,args)

