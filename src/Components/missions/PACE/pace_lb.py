#!/usr/bin/env python
"""
Wrapper to call leo_sampler.py
"""

import os
import argparse
from MAPL            import Config

#------------------------------------ M A I N ------------------------------------

if __name__ == "__main__":
    # Defaults
    parser = argparse.ArgumentParser()

    parser.add_argument("prep_config",
                        help="prep config filename")

    parser.add_argument("-v", "--verbose",action="store_true",
                        help="Verbose mode (default=False).")    


    args = parser.parse_args()

    bin = "./leo_sampler.py"

    # Parse prep config
    # -----------------
    cf = Config(args.prep_config,delim=' = ')

    iso_t1   = cf('ISO_T1')
    iso_t2   = cf('ISO_T2')
    algo     = cf('ALGO')
    rcFiles  = cf('RCFILES')

    algo = algo.split(',')
    rcFiles = rcFiles.split(',')

    for a,rc in zip(algo,rcFiles):
        # always write out coordinates
        command = bin + ' -C'

        if args.verbose:
            command += ' -v'

        command += ' --rcFile {}'.format(rc)
        command += ' --algo {}'.format(a)
        command += ' {} {}'.format(iso_t1,iso_t2)

        print command
        os.system(command)