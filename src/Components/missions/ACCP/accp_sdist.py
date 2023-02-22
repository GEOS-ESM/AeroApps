#!/usr/bin/env python3

"""
    Calculates the column size distribution for
    the GEOS sampled lidar track
"""
import os
import argparse
from   datetime        import datetime, timedelta
from   dateutil.parser import parse         as isoparser
from   MAPL            import Config
from   netCDF4 import Dataset
import numpy   as np
from   mieobs  import getEdgeVars
from MAPL.constants import *
from sdist import SDIST

# Generic Lists of Varnames and Units
VNAMES_DU = ['DU001','DU002','DU003','DU004','DU005']
VNAMES_SS = ['SS001','SS002','SS003','SS004','SS005']
VNAMES_BC = ['BCPHOBIC','BCPHILIC']
VNAMES_OC = ['OCPHOBIC','OCPHILIC']
VNAMES_SU = ['SO4']

META    = ['DELP','PS','RH','AIRDENS','LONGITUDE','LATITUDE','isotime']
AERNAMES = VNAMES_SU + VNAMES_SS + VNAMES_OC + VNAMES_BC + VNAMES_DU
AERdistNAMES = ['DU','SS','OCPHILIC','OCPHOBIC','BCPHILIC','BCPHOBIC','SU']
SDS_AER = META + AERNAMES

ncALIAS = {'LONGITUDE': 'trjLon',
           'LATITUDE' : 'trjLat'}

MISSING = -1.e+20

class ACCP_SDIST(SDIST):
    """
    Everything you need to create files of column size distribution
    """
    def __init__(self,inFile,outFile,rcFile,verbose=False):
        self.SDS_AER = SDS_AER
        self.AERNAMES = AERNAMES
        self.AERdistNAMES = AERdistNAMES
        self.inFile   = inFile
        self.outFile  = outFile
        self.rcFile   = rcFile
        self.verbose  = verbose

        # initialize empty lists
        for sds in self.SDS_AER:
            self.__dict__[sds] = []

        # Read in model data
        self.readSampledGEOS()

        # Make lists into arrays
        for sds in self.SDS_AER:
            self.__dict__[sds] = np.concatenate(self.__dict__[sds])

        # Get ze variable
        pe, ze, te = getEdgeVars(self)

        self.pe = pe # (km,nobs)
        self.ze = ze
        self.te = te


        # convert isotime to datetime
        self.tyme = []
        for isotime in self.isotime:
            self.tyme.append(isoparser(''.join(isotime)))

        self.tyme = np.array(self.tyme)
        self.ntyme  = len(self.tyme)

        # calculate column size distribution
        self.sizeDistribution()
 
        # write outfile
        self.writenc()
    #---
    def readSampledGEOS(self):
        """
        Read in model sampled track
        """
        col = 'aer_Nv'
        if self.verbose:
            print('opening file',self.inFile.replace('%col',col))
        nc       = Dataset(self.inFile.replace('%col',col))

        for sds in self.SDS_AER:
            sds_ = sds
            if sds in ncALIAS:
                sds_ = ncALIAS[sds]
            # make this explicitly array because 
            # in sdist.py interpolate doesn't handle masked arrays
            var = np.array(nc.variables[sds_][:])
            self.__dict__[sds].append(var)


if __name__ == "__main__":

    # Defaults
    DT_hours   = 1
    rcFile     = 'Aod_EOS.rc'

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

    parser.add_argument("--rcfile",default=rcFile,
                        help="rcFile (default=%s)"%rcFile)

    parser.add_argument("-D","--DT_hours", default=DT_hours, type=int,
                        help="Timestep in hours for each file (default=%i)"%DT_hours)

    parser.add_argument("-v", "--verbose",action="store_true",
                        help="Verbose mode (default=False).")

    parser.add_argument("-r", "--dryrun",action="store_true",
                        help="do a dry run (default=False).")

    args = parser.parse_args()
    rcFile         = args.rcfile

    # Parse prep config
    # -----------------
    cf             = Config(args.orbit_pcf,delim=' = ')
    orbitname      = cf('orbitname')
    ORBITNAME      = orbitname.upper()

    cf             = Config(args.track_pcf,delim=' = ')
    inTemplate     = cf('inDir')     + '/' + cf('inFile')
    outTemplate    = inTemplate

    # Loop through dates, calculating size distribution
    # -------------------------------------------------
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
        outFile    = outTemplate.replace('%year',year).replace('%month',month).replace('%day',day).replace('%nymd',nymd).replace('%hour',hour).replace('%orbitname',orbitname).replace('%ORBITNAME',ORBITNAME).replace('%col','aer_SD')

        # Initialize VLIDORT class getting aerosol optical properties
        # -----------------------------------------------------------
        print('++++Running VLIDORT with the following arguments+++')
        print('>>>inFile:    ',inFile)
        print('>>>outFile:   ',outFile)
        print('>>>rcFile:    ',rcFile)
        print('>>>verbose:   ',args.verbose)
        print('++++End of arguments+++')
        if not args.dryrun:
            sdist = ACCP_SDIST(inFile,outFile,rcFile,verbose=args.verbose)
            sdist = None


        date += Dt
