#!/usr/bin/env python
"""
Wrapper to call run_g5nr_stn_sampler.py for myd13c2
"""

import os
import argparse

#------------------------------------ M A I N ------------------------------------

if __name__ == "__main__":
    # Defaults
    DT_hours = 1

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


    args = parser.parse_args()

    command = "run_g5nr_stn_sampler.py"

    command += ' -D {}'.format(args.DT_hours)
    if args.verbose:
        command += ' -v'

    command += ' -a nearest {} {} {}'.format(args.iso_t1,args.iso_t2,args.prep_config)

    print command
    os.system(command)