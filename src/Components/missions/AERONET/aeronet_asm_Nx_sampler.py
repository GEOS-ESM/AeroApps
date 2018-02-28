#!/usr/bin/env python
"""
Wrapper to call run_g5nr_stn_sampler_invariant.py for asm_Nx
"""

import os
import argparse

#------------------------------------ M A I N ------------------------------------

if __name__ == "__main__":
    # Defaults
    parser = argparse.ArgumentParser()

    parser.add_argument("prep_config",
                        help="prep config filename")

    parser.add_argument("-v", "--verbose",action="store_true",
                        help="Verbose mode (default=False).")    


    args = parser.parse_args()

    command = "run_g5nr_stn_sampler_invariant.py"

    if args.verbose:
        command += ' -v'

    command += ' {}'.format(args.prep_config)

    print command
    os.system(command)