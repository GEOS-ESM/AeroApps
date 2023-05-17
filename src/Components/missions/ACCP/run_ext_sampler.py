#!/usr/bin/env python3

"""
    Wrapper for ext_sampler.py
    Gets extinction veritcal profiles and aerosol intensive properties
    Uses trj_sampled files as inputs

    adapted from run_lidar_sampler.py
    Patricia Castellanos, Oct, 2019

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
      print(cmds[nextdate])
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
    nproc    = 12
    rcFile   = 'Aod_EOS.rc'

    parser = argparse.ArgumentParser()
    parser.add_argument("iso_t1",
                        help="starting iso time")
    parser.add_argument("iso_t2",
                        help="ending iso time")

    parser.add_argument("input",
                        help="filename template")

    parser.add_argument("output",
                        help="filename template")   

    parser.add_argument("channel",
                        help="channel")                            

    parser.add_argument("-D","--DT_hours", default=DT_hours, type=int,
                        help="Timestep in hours for each file (default=%i)"%DT_hours)

    parser.add_argument("--rc", default=rcFile,
                        help="rc File (default=%s)"%rcFile)    

    parser.add_argument("-v", "--verbose",action="store_true",
                        help="Verbose mode (default=False).")

    parser.add_argument("-r", "--dryrun",action="store_true",
                        help="do a dry run (default=False).")    

    parser.add_argument("--coarse",action="store_true",
                        help="do coarse mode only (default=False).")

    parser.add_argument("--fine",action="store_true",
                        help="do fine mode only (default=False).")
    parser.add_argument("--spc",action="store_true",
                        help="do indivudual species (default=False).")

    parser.add_argument("-n", "--nproc",default=nproc,type=int,
                        help="Number of processors (default=%i)."%nproc)    

    args = parser.parse_args()

    ch        = args.channel   
    enddate   = isoparser(args.iso_t2)
    pdf = timedelta(hours=DT_hours)

    Date = isoparser(args.iso_t1)
    inFilelist = []
    outFilelist = []
    while Date < enddate:
        YY = Date.strftime('%Y')
        MM = Date.strftime('%m')
        DD = Date.strftime('%d')
        strdate = Date.strftime('%Y%m%d_%H%Mz')            
        inFile = args.input.format(YY,MM,DD,strdate)
        outFile = args.output.format(YY,MM,DD,strdate,ch)
        outpath = os.path.dirname(outFile)
        mkdir_p(outpath)

        inFilelist.append(inFile)
        outFilelist.append(outFile)

        Date += timedelta(hours=args.DT_hours)

    # run extinction sampler on model fields
    # split across multiple processors by date
    processes = []
    cmds      = []
    for inFile,outFile in zip(inFilelist,outFilelist):
        # all species
        Options =     " --vnames=SS001,SS002,SS003,SS004,SS005,BCPHOBIC,BCPHILIC,OCPHOBIC,OCPHILIC,SO4,DU001,DU002,DU003,DU004,DU005"     + \
                      " --input=" + inFile       + \
                      " --output=" + outFile      + \
                      " --channel=" + ch      + \
                      " --rc=" + args.rc      + \
                      " --format=NETCDF4_CLASSIC"      + \
                      " --intensive"  

        if args.verbose:
            Options += " --verbose" 

        cmd = 'ext_sampler.py {} '.format(Options)
        cmds.append(cmd)

    if args.spc:
        SPC = {'SU':'SO4',
               'OC':'OCPHILIC,OCPHOBIC',
               'BC':'BCPHILIC,BCPHOBIC',
               'SS':'SS001,SS002,SS003,SS004,SS005',
               'DU':'DU001,DU002,DU003,DU004,DU005'}
        for spc in SPC:
            for inFile,outFile in zip(inFilelist,outFilelist):
            # all species indivudually    
                Options = " --vnames=" + SPC[spc] + \
                              " --input=" + inFile       + \
                              " --output=" + outFile[:-4] + "."+spc+".nc4"      + \
                              " --channel=" + ch      + \
                              " --rc=" + args.rc      + \
                              " --format=NETCDF4_CLASSIC" + \
                              " --intensive"

                if args.verbose:
                    Options += " --verbose"

                cmd = 'ext_sampler.py {} '.format(Options)
                cmds.append(cmd)

    # coarse mode only
    if args.coarse:
        for inFile,outFile in zip(inFilelist,outFilelist):
            Options =     " --vnames=SS001,SS002,SS003,SS004,SS005,DU001,DU002,DU003,DU004,DU005"     + \
                              " --input=" + inFile       + \
                              " --output=" + outFile[:-4] + ".coarse.nc4"      + \
                              " --channel=" + ch      + \
                              " --rc=" + args.rc      + \
                              " --format=NETCDF4_CLASSIC"      + \
                              " --intensive"

            if args.verbose:
                Options += " --verbose"

            cmd = 'ext_sampler.py {} '.format(Options)
            cmds.append(cmd)

    # fine mode only
    if args.fine:
        for inFile,outFile in zip(inFilelist,outFilelist):
            Options =     " --vnames=BCPHOBIC,BCPHILIC,OCPHOBIC,OCPHILIC,SO4"     + \
                              " --input=" + inFile       + \
                              " --output=" + outFile[:-4] + ".fine.nc4"      + \
                              " --channel=" + ch      + \
                              " --rc=" + args.rc      + \
                              " --format=NETCDF4_CLASSIC"      + \
                              " --intensive"

            if args.verbose:
                Options += " --verbose"

            cmd = 'ext_sampler.py {} '.format(Options)
            cmds.append(cmd)

    lendate  = len(cmds)
    # Manage processes
    # This will start the max processes running    
    if not args.dryrun:
        processes, nextdate = CheckRunning(processes,cmds,0,lendate,args)
        while len(processes)>0: # Some things still going on
            time.sleep(10)      # Wait
            # add more processes as other ones finish
            processes, nextdate = CheckRunning(processes,cmds,nextdate,lendate,args)
    else:
        for cmd in cmds:
            print(cmd)

