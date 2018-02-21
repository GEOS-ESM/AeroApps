#!/usr/bin/env python

"""
    Wrapper for run_aeronet_vlidort.py
    Runs vlidort simulator

    Patricia Castellanos, May, 2017

"""

import os
import argparse
from   datetime        import datetime, timedelta
from   dateutil.parser import parse         as isoparser
import numpy  as np

#------------------------------------ M A I N ------------------------------------

if __name__ == "__main__":

    # Defaults
    DT_hours = 24
    channels = ['340','380','440','500','670','870','1020','1640']
    extch    = ['355','532','1064']
    platform = 'nccs'

#   Parse command line options
#   --------------------------

    parser = argparse.ArgumentParser()
    parser.add_argument("iso_t1",
                        help="starting iso time")
    parser.add_argument("iso_t2",
                        help="ending iso time")

    parser.add_argument("-D","--DT_hours", default=DT_hours, type=int,
                        help="Timestep in hours for each file (default=%i)"%DT_hours)

    parser.add_argument("-v", "--verbose",action="store_true",
                        help="Verbose mode (default=False).")

    parser.add_argument("-r", "--dryrun",action="store_true",
                        help="do a dry run (default=False).")   

    parser.add_argument("-c","--channels", default=channels,nargs='+',
                        help="channels to get BOA radiance")                            

    parser.add_argument("-e","--extch", default=extch,nargs='+',
                        help="channels to run extinction sampler")  

    parser.add_argument("-p","--platform", default=platform,
                        help="options are 'calculon' or 'nccs'")  

    parser.add_argument("--norad",action="store_true",
                        help="No radiance calculations (default=False).")

    parser.add_argument("--noext",action="store_true",
                        help="No extinctions calculations (default=False).")


    args = parser.parse_args()


    # Loop through dates and channels
    # ------------------------------------
    bin = './run_aeronet_vlidort.py'

    if not args.norad:
        for ch in args.channels:
            # Run aeronet_vlidort
            # -----------------------------------------------------------
            print '++++Running VLIDORT {} nm+++'.format(ch)

            if args.dryrun:
                command = bin + '-r '
            else:
                command = bin


            command = command + ' -D {} {} {} {}/stn_vlidort_{}.pcf'.format(args.DT_hours,args.iso_t1,args.iso_t2,args.platform,ch)


            print command
            os.system(command)

    if not args.noext:
        for ch in args.extch:
            # Run extinction sampler
            # -----------------------------------------------------------
            print '++++Running EXT_SAMPLER {} nm+++'.format(ch)

            if args.dryrun:
                command = bin + '-r '
            else:
                command = bin


            command = command + ' -D {} {} {} {}/stn_vlidort_{}ext.pcf'.format(args.DT_hours,args.iso_t1,args.iso_t2,args.platform,ch)


            print command
            os.system(command)


