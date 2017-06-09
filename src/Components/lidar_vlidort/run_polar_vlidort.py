#!/usr/bin/env python

"""
    Wrapper for polar_vlidort.py
    Runs vlidort simulator

    Patricia Castellanos, May, 2017

"""

import os
import argparse
from   datetime        import datetime, timedelta
from   dateutil.parser import parse         as isoparser
from   MAPL            import Config
from   polar_vlidort   import POLAR_VLIDORT, get_chd
import numpy  as np

#------------------------------------ M A I N ------------------------------------

if __name__ == "__main__":

    # Defaults
    DT_hours = 24

#   Parse command line options
#   --------------------------

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

    args = parser.parse_args()

    # Parse prep config
    # -----------------
    cf             = Config(args.prep_config,delim=' = ')
    inTemplate     = cf('inDir')     + '/' + cf('inFile')        
    albedoTemplate = cf('albedoDir') + '/' + cf('albedoFile')  
    outTemplate    = cf('outDir')    + '/' + cf('outFile')
    channel        = int(cf('channel'))
    rcFile         = cf('rcFile')

    VZA           = cf('VZA')
    VZA = np.array(VZA.replace(' ','').split(',')).astype('float')

    if cf('ndviDir').upper() == 'NONE':
        ndviTemplate = None
        lcTemplate   = None
    else:
        ndviTemplate   = cf('ndviDir')   + '/' + cf('ndviFile')
        lcTemplate     = cf('lcDir')     + '/' + cf('lcFile')



    # Change wavelength number in Aod_EOS.rc
    # ---------------------------------------
    source = open(rcFile,'r')
    destination = open(rcFile+'.'+cf('channel')+'.tmp','w')
    a = float(channel)*1e-3
    for line in source:
        if (line[0:11] == 'r_channels:'):
            destination.write('r_channels: '+'{:0.3f}e-6'.format(a)+'\n')
        else:
            destination.write(line) 
    source.close()
    destination.close()
    rcFile = rcFile+'.'+cf('channel')+'.tmp'

    # Loop through dates, running VLIDORT
    # ------------------------------------
    date      = isoparser(args.iso_t1)
    enddate   = isoparser(args.iso_t2)
    Dt        = timedelta(hours=args.DT_hours)

    while date < enddate:
        nymd = str(date.date()).replace('-','')
        year = str(date.year)
        month = str(date.month).zfill(2)
        hour = str(date.hour).zfill(2)    

        inFile     = inTemplate.replace('%year',year).replace('%month',month).replace('%nymd',nymd).replace('%hour',hour)
        outFile    = outTemplate.replace('%year',year).replace('%month',month).replace('%nymd',nymd).replace('%hour',hour)
        albedoFile = albedoTemplate.replace('%year',year).replace('%month',month).replace('%nymd',nymd).replace('%hour',hour)
        
        if ndviTemplate is None:
            ndviFile = None
            lcFile   = None
        else:
            ndviFile   = ndviTemplate.replace('%year',year).replace('%month',month).replace('%nymd',nymd).replace('%hour',hour)
            lcFile     = lcTemplate.replace('%year',year).replace('%month',month).replace('%nymd',nymd).replace('%hour',hour)

        # Initialize VLIDORT class getting aerosol optical properties
        # -----------------------------------------------------------
        print '++++Running VLIDORT with the following arguments+++'
        print 'inFile:',inFile
        print 'outFile:',outFile
        print 'rcFile:',rcFile
        print 'albedoFile:',albedoFile
        print 'albedoType:',cf('albedoType')
        print 'channel:',channel
        print 'VZA:',VZA
        print 'HGT:',float(cf('HGT'))
        print 'ndviFile:',ndviFile,
        print 'lcFile:',lcFile,
        print 'verbose:',args.verbose
        if not args.dryrun:
            vlidort = POLAR_VLIDORT(inFile,outFile,rcFile,
                                    albedoFile,cf('albedoType'),
                                    channel,
                                    VZA,
                                    float(cf('HGT')),
                                    ndviFile=ndviFile,
                                    lcFile=lcFile,
                                    verbose=args.verbose)

            # Run VLIDORT
            if vlidort.nobs > 0:
                if cf('DO_VLIDORT').upper() == 'YES':
                    vlidort.runVLIDORT()
                if cf('DO_EXT').upper() == 'YES':
                    vlidort.runExt()

        date += Dt


    # clean up rcFile
    os.remove(rcFile)
