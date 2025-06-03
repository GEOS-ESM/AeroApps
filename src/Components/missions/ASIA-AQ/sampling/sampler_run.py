#!/usr/bin/env python3
"""
    Driver script for aaq_sampler.py
"""
import os, sys
import subprocess
import argparse
import yaml
from glob import glob
import time

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

   while (len(processes) < int(args.nproc)) and (nextdate < lendate): # More to do and some spare slots
      processes, nextdate = StartNew(processes,cmds,nextdate,lendate)

   return processes,nextdate


if __name__ == "__main__":


    parser = argparse.ArgumentParser()
    parser.add_argument("config",help='configuration yaml file')
    parser.add_argument("--nproc",default=1)

    args = parser.parse_args()

    config = yaml.safe_load(open(args.config))

    # create output directory
    if not os.path.exists(config['sampled_outdir']):
        os.makedirs(config['sampled_outdir'])

    # get dc-8 files
    ictFiles = sorted(glob(config['dc8_merge']+'/*ict'))

    cmds = []
    for ict in ictFiles:
        cmd = f'./aaq_sampler.py {ict} {args.config}'
        cmds.append(cmd)
    
    lendate = len(cmds)    
    # Manage processes
    # This will start the max processes running
    processes = []
    processes, nextdate = CheckRunning(processes,cmds,0,lendate,args)
    while len(processes)>0: # Some things still going on
        time.sleep(10)      # Wait
        # add more processes as other ones finish
        processes, nextdate = CheckRunning(processes,cmds,nextdate,lendate,args)

