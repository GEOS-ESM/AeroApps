#!/usr/bin/env python3

"""
    Wrapper for lidar_vlidort.py
    Runs vlidort simulator

    Patricia Castellanos, Dec, 2019

"""

import os
import argparse
from   datetime        import datetime, timedelta
from   dateutil.parser import parse         as isoparser
from   MAPL.config     import Config
from   py_leo_vlidort.lidar_vlidort   import LIDAR_VLIDORT, get_chd
import numpy  as np

#------------------------------------ M A I N ------------------------------------

if __name__ == "__main__":

    # Defaults
    DT_hours   = 1
    instname   = 'lidar'
    rcFile     = 'Aod_EOS.rc'
    albedoType = None

#   Parse command line options
#   --------------------------

    parser = argparse.ArgumentParser()
    parser.add_argument("iso_t1",
                        help="starting iso time")
    parser.add_argument("iso_t2",
                        help="ending iso time")

    parser.add_argument("track_pcf",
                        help="prep config file with track input file names")

    parser.add_argument("orbit_pcf",
                        help="prep config file with orbit variables")

    parser.add_argument("channel", type=int,
                        help="channel in nm")

    parser.add_argument("-a","--albedotype", default=albedoType,
                        help="albedo type keyword. default is to figure out according to channel")

    parser.add_argument("-I","--instname", default=instname,
                        help="Instrument name (default=%s)"%instname)
    parser.add_argument("--rcfile",default=rcFile,
                        help="rcFile (default=%s)"%rcFile)

    parser.add_argument("-D","--DT_hours", default=DT_hours, type=int,
                        help="Timestep in hours for each file (default=%i)"%DT_hours)

    parser.add_argument("-v", "--verbose",action="store_true",
                        help="Verbose mode (default=False).")

    parser.add_argument("-r", "--dryrun",action="store_true",
                        help="do a dry run (default=False).")    

    args = parser.parse_args()
    channel        = args.channel
    rcFile         = args.rcfile
    instname       = args.instname
    albedoType     = args.albedotype

    # figure out albedoType keyword
    if albedoType is None:
        if channel <= 388:
            albedoType = 'LAMBERTIAN'
        else:
            albedoType = 'MODIS_BRDF'

    # Parse prep config
    # -----------------
    cf             = Config(args.orbit_pcf,delim=' = ')
    orbitname      = cf('orbitname')
    ORBITNAME      = orbitname.upper()
    HGT            = float(cf('HGT'))

    cf             = Config(args.track_pcf,delim=' = ')
    inTemplate     = cf('inDir')     + '/' + cf('inFile')         
    outTemplate    = cf('outDir')    + '/' + cf('outFile')
    
    try:
        brdfTemplate = cf('brdfDir') + '/' + cf('brdfFile') 
    except:
        brdfTemplate = None

    try:
        ndviTemplate   = cf('ndviDir')   + '/' + cf('ndviFile')
        lcTemplate     = cf('lcDir')     + '/' + cf('lcFile')
    except:
        ndviTemplate = None
        lcTemplate   = None

    try:
        lerTemplate    = cf('lerDir')    + '/' + cf('lerFile')
    except:
        lerTemplate    = None


    # Change wavelength number in Aod_EOS.rc
    # ---------------------------------------
    source = open(rcFile,'r')
    destination = open(rcFile+'.'+str(channel)+'.tmp','w')
    a = float(channel)*1e-3
    for line in source:
        if (line[0:11] == 'r_channels:'):
            destination.write('r_channels: '+'{:0.3f}e-6'.format(a)+'\n')
        else:
            destination.write(line) 
    source.close()
    destination.close()
    rcFile = rcFile+'.'+str(channel)+'.tmp'

    # Loop through dates, running VLIDORT
    # ------------------------------------
    date      = isoparser(args.iso_t1)
    enddate   = isoparser(args.iso_t2)
    Dt        = timedelta(hours=args.DT_hours)

    while date < enddate:
        nymd  = str(date.date()).replace('-','')
        year  = str(date.year)
        month = str(date.month).zfill(2)
        day   = str(date.day).zfill(2)
        hour  = str(date.hour).zfill(2)    

        inFile     = inTemplate.replace('%year',year).replace('%month',month).replace('%day',day).replace('%nymd',nymd).replace('%hour',hour).replace('%orbitname',orbitname).replace('%ORBITNAME',ORBITNAME)
        outFile    = outTemplate.replace('%year',year).replace('%month',month).replace('%day',day).replace('%nymd',nymd).replace('%hour',hour).replace('%chd',get_chd(channel)).replace('%orbitname',orbitname).replace('%ORBITNAME',ORBITNAME).replace('%instname',instname)

        if brdfTemplate is None:
            brdfFile = None
        else:
            brdfFile = brdfTemplate.replace('%year',year).replace('%month',month).replace('%day',day).replace('%nymd',nymd).replace('%hour',hour).replace('%orbitname',orbitname).replace('%ORBITNAME',ORBITNAME)
        
        if ndviTemplate is None:
            ndviFile = None
            lcFile   = None
        else:
            ndviFile   = ndviTemplate.replace('%year',year).replace('%month',month).replace('%day',day).replace('%nymd',nymd).replace('%hour',hour).replace('%orbitname',orbitname).replace('%ORBITNAME',ORBITNAME)
            lcFile     = lcTemplate.replace('%year',year).replace('%month',month).replace('%day',day).replace('%nymd',nymd).replace('%hour',hour).replace('%orbitname',orbitname).replace('%ORBITNAME',ORBITNAME)

        if lerTemplate is None:
            lerFile = None
        else:
            lerFile   = lerTemplate.replace('%year',year).replace('%month',month).replace('%day',day).replace('%nymd',nymd).replace('%hour',hour).replace('%orbitname',orbitname).replace('%ORBITNAME',ORBITNAME)

        # Initialize VLIDORT class getting aerosol optical properties
        # -----------------------------------------------------------
        print('++++Running VLIDORT with the following arguments+++')
        print('>>>inFile:    ',inFile)
        print('>>>outFile:   ',outFile)
        print('>>>rcFile:    ',rcFile)
        print('>>>albedoType:',albedoType)
        print('>>>channel:   ',channel)
        print('>>>HGT:       ',HGT)
        print('>>>brdfFile:  ',brdfFile)
        print('>>>ndviFile:  ',ndviFile)
        print('>>>lcFile:    ',lcFile)
        print('>>>lerFile    ',lerFile)
        print('>>>verbose:   ',args.verbose)
        print('++++End of arguments+++')
        if not args.dryrun:
            vlidort = LIDAR_VLIDORT(inFile,outFile,rcFile,
                                    albedoType,
                                    channel,
                                    HGT,
                                    brdfFile=brdfFile,
                                    ndviFile=ndviFile,
                                    lcFile=lcFile,
                                    lerFile=lerFile,
                                    verbose=args.verbose)

            # Run VLIDORT
            if vlidort.nobs > 0:
                vlidort.runVLIDORT()

        date += Dt


    # clean up rcFile
    os.remove(rcFile)
