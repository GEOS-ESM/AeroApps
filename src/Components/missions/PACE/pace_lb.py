#!/usr/bin/env python
"""
Wrapper to do all sampling leo_sampler.py
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

    # Parse prep config
    # -------------------
    cf = Config(args.prep_config,delim=' = ')

    iso_t1   = cf('ISO_T1')
    iso_t2   = cf('ISO_T2')


    # G5NR Sampling
    # -------------------
    try:
        g5nrPCF = cf('G5NR')

    except:
        g5nrPCF = None

    if g5nrPCF is not None:
        bin = "./leo_sampler.py"

        # Parse prep config
        # -----------------
        cf = Config(g5nrPCF,delim=' = ')

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

    # MCD12C (Land Cover Type) Sampling 
    # -------------------------------------
    try:
        mcd12cPCF = cf('MCD12C')

    except:
        MCD12CPCF = None   

    if mcd12cPCF is not None:
        bin = "./mcd12c_sampler.py"

        command = ''
        if args.verbose:
            command += ' -v'

        command += ' {} {} {}'.format(iso_t1,iso_t2,mcd12cPCF)

        print command
        os.system(command)

    # MCD43C (BRDF) Sampling
    # --------------------------
    try:
        mcd43cPCF = cf('MCD43C')

    except:
        mcd43cPCF = None   

    if mcd43cPCF is not None:
        bin = "./mcd43c_sampler.py"

        command = ''
        if args.verbose:
            command += ' -v'

        command += ' {} {} {}'.format(iso_t1,iso_t2,mcd43cPCF)

        print command
        os.system(command)

    # NOBM (Water Leaving Radiance) Sampling
    # --------------------------
    try:
        nobmPCF = cf('NOBM')

    except:
        nobmPCF = None   

    if nobmPCF is not None:
        bin = "./nobm_sampler.py"

        command = ''
        if args.verbose:
            command += ' -v'

        command += ' {} {} {}'.format(iso_t1,iso_t2,nobmPCF)

        print command
        os.system(command)


    # MYD13C2 (NDVI) Sampling
    # ------------------------
    try:
        myd13c2PCF = cf('MYD13C2')
    except:
        myd13c2PCF = None

    if myd13c2PCF is not None:
        bin = "./leo_sampler.py"

        # always write out coordinates
        command = bin + ' -C'

        if args.verbose:
            command += ' -v'

        command += ' --rcFile {}'.format(myd13c2PCF)
        command += ' --algo nearest'
        command += ' {} {}'.format(iso_t1,iso_t2)

        print command
        os.system(command)

    # OMI LER Sampling
    # ----------------------
    try:
        lerPCF = cf('LER')
    except:
        lerPCF = None

    if myd13c2PCF is not None:
        bin = "./leo_sampler.py"

        # always write out coordinates
        command = bin + ' -C'

        if args.verbose:
            command += ' -v'

        command += ' --rcFile {}'.format(lerPCF)
        command += ' --algo linear'
        command += ' {} {}'.format(iso_t1,iso_t2)

        print command
        os.system(command)




