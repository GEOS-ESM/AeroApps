#!/usr/bin/env python

"""
    Puts all needed lidar signal simulator variables in one file

    Patricia Castellanos, Dec, 2019

"""

import os
import argparse
from   datetime        import datetime, timedelta
from   dateutil.parser import parse         as isoparser
from   MAPL            import Config
import numpy  as np
from netCDF4 import Dataset
import LidarAngles_   

# Generic Lists of Varnames and Units
SDS_AER    = ['LONGITUDE','LATITUDE','isotime']


ncALIAS = {'LONGITUDE': 'trjLon',
           'LATITUDE': 'trjLat'}


class SWATH(object):
    """
    get polarimeter swath
    """
    def __init__(self,inFile,outFile,hgtss):
        self.inFile = inFile
        self.outFile = outFile
        self.hgtss     = hgtss
        self.SDS_AER    = SDS_AER

        # initialize empty lists
        for sds in self.SDS_AER:
            self.__dict__[sds] = []

        # Read in model data
        self.readSampledGEOS()

        # Make lists into arrays
        for sds in self.SDS_AER:
            self.__dict__[sds] = np.concatenate(self.__dict__[sds])

        # convert isotime to datetime
        self.tyme = []
        for isotime in self.isotime:
            self.tyme.append(isoparser(''.join(isotime)))

        self.tyme = np.array(self.tyme)

        # Calculate Nadir Scene Geometry
        self.hgtss = hgtss
        self.nadirAngles()
    # --
    def nadirAngles(self):
        SZA   = []
        SAA   = []
        VAA   = []

        for i,tyme in enumerate(self.tyme[:-2]):
            CLAT = self.LATITUDE[i]
            CLON = self.LONGITUDE[i]
            
            SLAT = self.LATITUDE[i]
            SLON = self.LONGITUDE[i]
            sat_angles = LidarAngles_.satangles(tyme.year,tyme.month,tyme.day,
                                                tyme.hour,tyme.minute,tyme.second,
                                                CLAT,CLON,
                                                SLAT,SLON,
                                                0.0,
                                                self.hgtss)
            SZA.append(sat_angles[3][0])
            SAA.append(sat_angles[2][0])


            CLAT = self.LATITUDE[i+1]
            CLON = self.LONGITUDE[i+1]
            sat_angles = LidarAngles_.satangles(tyme.year,tyme.month,tyme.day,
                                                tyme.hour,tyme.minute,tyme.second,
                                                CLAT,CLON,
                                                SLAT,SLON,
                                                0.0,
                                                self.hgtss)

            VAA.append(sat_angles[0][0])


        self.SZA = np.array(SZA)
        self.SAA = np.array(SAA)
        self.VAA = np.array(VAA)     

    #---
    def readSampledGEOS(self):
        """
        Read in model sampled track
        """
        col = 'aer_Nv'
        if self.verbose: 
            print 'opening file',self.inFile.replace('%col',col)
        nc       = Dataset(self.inFile.replace('%col',col))

        for sds in self.SDS_AER:
            sds_ = sds
            if sds in ncALIAS:
                sds_ = ncALIAS[sds]
            var = nc.variables[sds_][:]
            self.__dict__[sds].append(var)



#------------------------------------ M A I N ------------------------------------

if __name__ == "__main__":

    # Defaults
    DT_hours = 1

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
    instname       = cf('instname')
    INSTNAME       = instname.upper()
    HGT            = float(cf('HGT'))
    polarimeter_name       = cf('polarimeter_name')


    # Loop through dates, calculating geometry
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

        inFile     = inTemplate.replace('%year',year).replace('%month',month).replace('%day',day).replace('%nymd',nymd).replace('%hour',hour).replace('%instname',instname).replace('%INSTNAME',INSTNAME)
        outFile    = inTemplate.replace('%col',polarimeter_name).replace('%year',year).replace('%month',month).replace('%day',day).replace('%nymd',nymd).replace('%hour',hour).replace('%instname',instname).replace('%INSTNAME',INSTNAME)


        # Initialize MERGE class and write new outfile
        # -----------------------------------------------------------
        print '++++Getting polarimeter swath with the following arguments+++'
        print '>>>inFile:    ',inFile
        print '>>>outFile:   ',outFile
        print '>>>HGT:       ',HGT
        print '>>>verbose:   ',args.verbose
        print '++++End of arguments+++'
        if not args.dryrun:
            swath = SWATH(inFile,outFile,HGT)