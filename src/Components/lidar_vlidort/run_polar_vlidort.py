#!/usr/bin/env python

"""
    Wrapper for polar_vlidort.py
    Runs vlidort simulator

    Patricia Castellanos, May, 2017

"""

import os
import argparse
from datetime        import datetime, timedelta
from dateutil.parser import parse         as isoparser
from MAPL            import Config
from polar_vlidort   import POLAR_VLIDORT

#------------------------------------ M A I N ------------------------------------

if __name__ == "__main__":

    # Defaults
    DT_hours = 24

#   Parse command line options
#   --------------------------
    parser = OptionParser()

    parser = argparse.ArgumentParser()
    parser.add_argument("iso_t1",
                        help="starting iso time")
    parser.add_argument("iso_t2",
                        help="ending iso time")

    parser.add_argument("prep_config",
                        help="prep config filename")

    parser.add_argument("-D","--DT_hours", default=DT_hours, type=int,
                        help="Timestep in hours for each file (default=%i)"%DT_hours)

    parser.add_argument("-v", "--verbose",action="store_true",
                        help="Verbose mode (default=False).")

    parser.add_argument("-r", "--dryrun",action="store_true",
                        help="do a dry run (default=False).")    

    parser.add("-c", "--channel", dest="channel", default=channel,type="int",
              help="Channel for Mie calculation")

    parser.add_option("-v", "--verbose",
                      action="store_true", dest="verbose",
                      help="Verbose mode")

    (options, args) = parser.parse_args()


    # edit Aod_Eos.rc

         
    # Create consistent file name extension
    # -------------------------------------
    name, ext = os.path.splitext(options.outFile)
    if 'NETCDF4' in options.format:
        options.outFile = name + '.nc4'
    elif 'NETCDF3' in options.format:
        options.outFile = name + '.nc'
    else:
        raise ValueError, 'invalid extension <%s>'%ext
