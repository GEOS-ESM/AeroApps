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



    args = parser.parse_args()


    # Loop through dates and channels
    # ------------------------------------
    date      = isoparser(args.iso_t1)
    enddate   = isoparser(args.iso_t2)
    Dt        = timedelta(hours=args.DT_hours)

    while date < enddate:
        for ch in args.channels:
            # Initialize VLIDORT class getting aerosol optical properties
            # -----------------------------------------------------------
            print '++++Running VLIDORT {} {} nm+++'.format(date.strftime('%Y%m%d %H'),ch)

            command = './run_aeronet_vlidort.py -D {} '.format(args.DT_hours)
            if args.dryrun:
                command = command + '-r '


            command = command + '{} {} {}/stn_vlidort_{}.pcf'.format(date.isoformat(),date.isoformat(),args.platform,ch)


            print command
            os.system(command)

        date += Dt