#!/usr/bin/env python3
"""
Wrapper to do all sampling leo_sampler.py
"""

import os
import argparse
from MAPL            import Config
from datetime        import datetime, timedelta
from dateutil.parser import parse         as isoparser

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

    piso_t1   = cf('ISO_T1')
    piso_t2   = cf('ISO_T2')

    t1 = isoparser(piso_t1)
    t2 = isoparser(piso_t2)

    giso_t1 = t1.strftime('2006-%m-%dT%H:%M:00')
    giso_t2 = t2.strftime('2006-%m-%dT%H:%M:00')

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
        cfg5nr = Config(g5nrPCF,delim=' = ')

        algo     = cfg5nr('ALGO')
        rcFiles  = cfg5nr('RCFILES')

        algo = algo.split(',')
        rcFiles = rcFiles.split(',')

        for a,rc in zip(algo,rcFiles):
            # always write out coordinates
            command = bin + ' -C'

            if args.verbose:
                command += ' -v'

            command += ' --rcFile {}'.format(rc)
            command += ' --algo {}'.format(a)
            command += ' {} {}'.format(piso_t1,piso_t2)

            print(command)
            os.system(command)

    # MCD12C (Land Cover Type) Sampling 
    # -------------------------------------
    try:
        mcd12cPCF = cf('MCD12C')

    except:
        mcd12cPCF = None   

    if mcd12cPCF is not None:
        bin = "./mcd12c_sampler.py"

        command = bin + ''
        if args.verbose:
            command += ' -v'

        command += ' {} {} {}'.format(giso_t1,giso_t2,mcd12cPCF)

        print(command)
        os.system(command)

    # MCD43C (BRDF) Sampling
    # --------------------------
    try:
        mcd43cPCF = cf('MCD43C')

    except:
        mcd43cPCF = None   

    if mcd43cPCF is not None:
        bin = "./mcd43c_sampler.py"

        command = bin + ''
        if args.verbose:
            command += ' -v'

        command += ' {} {} {}'.format(giso_t1,giso_t2,mcd43cPCF)

        print(command)
        os.system(command)

    # NOBM (Water Leaving Radiance) Sampling
    # --------------------------
    try:
        nobmPCF = cf('NOBM')

    except:
        nobmPCF = None   

    if nobmPCF is not None:
        bin = "./nobm_sampler.py"

        command = bin + ''
        if args.verbose:
            command += ' -v'

        command += ' {} {} {}'.format(giso_t1,giso_t2,nobmPCF)

        print(command)
        os.system(command)


    # MYD13C2 (NDVI) Sampling
    # ------------------------
    try:
        myd13c2PCF = cf('MYD13C2')
    except:
        myd13c2PCF = None

    if myd13c2PCF is not None:
        bin = "./myd13c2_sampler.py"

        command = bin + ''
        if args.verbose:
            command += ' -v'

        command += ' {} {} {}'.format(giso_t1,giso_t2,myd13c2PCF)

        print(command)
        os.system(command)

    # OMI LER Sampling
    # ----------------------
    try:
        lerPCF = cf('LER')
    except:
        lerPCF = None

    if myd13c2PCF is not None:
        bin = "./omi_ler_sampler.py"

        command = bin + ''
        if args.verbose:
            command += ' -v'

        command += ' {} {} {}'.format(giso_t1,giso_t2,lerPCF)

        print(command)
        os.system(command)




