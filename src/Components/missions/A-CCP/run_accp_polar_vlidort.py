#!/usr/bin/env python

"""
    Wrapper for accp_polar_vlidort.py
    Uses subprocess to run more than one instance

    Patricia Castellanos, Jan 2020

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

#------------------------------------ M A I N ------------------------------------
def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise

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
    DT_hours = 1
    nproc    = 8
    rcFile   = 'Aod_EOS.rc'
    albedoType = None

    parser = argparse.ArgumentParser()
    parser.add_argument("iso_t1",
                        help="starting iso time")
    parser.add_argument("iso_t2",
                        help="ending iso time")

    parser.add_argument("track_pcf",
                        help="prep config file with orbit variables")

    parser.add_argument("orbit_pcf",
                       help="prep config file with orbit variables")

    parser.add_argument("inst_pcf",
                        help="prep config file with instrument variables")

    parser.add_argument("channel", type=int,
                        help="channel in nm")

    parser.add_argument("-a","--albedotype", default=albedoType,
                        help="albedo type keyword. default is to figure out according to channel")

    parser.add_argument("--rcfile",default=rcFile,
                        help="rcFile (default=%s)"%rcFile)

    parser.add_argument("-D","--DT_hours", default=DT_hours, type=int,
                        help="Timestep in hours for each file (default=%i)"%DT_hours)

    parser.add_argument("-v", "--verbose",action="store_true",
                        help="Verbose mode (default=False).")

    parser.add_argument("-r", "--dryrun",action="store_true",
                        help="do a dry run (default=False).")    

    parser.add_argument("-n", "--nproc",default=nproc,type=int,
                        help="Number of processors (default=%i)."%nproc)    

    args = parser.parse_args()
    track_pcf      = args.track_pcf
    orbit_pcf      = args.orbit_pcf
    inst_pcf       = args.inst_pcf
    channel        = args.channel
    rcFile         = args.rcfile
    albedoType     = args.albedotype


    Date = isoparser(args.iso_t1)
    enddate   = isoparser(args.iso_t2)
    pdt   = timedelta(hours=args.DT_hours)
    DT    = enddate - Date
    DT    = int(DT.total_seconds()/(60*60))
    datelist = [Date + p*pdt for p in range(DT)]

    # run vlidort 
    # split across multiple processors by date
    lendate  = len(datelist)

    processes = []
    cmds      = []
    for date in datelist:
        edate = date + pdt 

        Options =     " --rcfile=" + rcFile      + \
                      " --DT_hours=" + str(args.DT_hours)

        if args.verbose:
            Options += " --verbose" 

        if args.dryrun:
            Options += " --dryrun"

        if albedoType is not None:
            Options += " --albedotype=" + albedoType

        cmd = './run_accp_polar_vlidort.py {} {} {} {} {} {} {}'.format(Options,date.isoformat(),edate.isoformat(),track_pcf,orbit_pcf,inst_pcf,channel)
        cmds.append(cmd)

    # Manage processes
    # This will start the max processes running    
    processes, nextdate = CheckRunning(processes,cmds,0,lendate,args)
    while len(processes)>0: # Some things still going on
        time.sleep(10)      # Wait
        # add more processes as other ones finish
        processes, nextdate = CheckRunning(processes,cmds,nextdate,lendate,args)


