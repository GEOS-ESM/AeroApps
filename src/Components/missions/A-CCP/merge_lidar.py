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
from   lidar_vlidort   import get_chd
import numpy  as np
from   netCDF4 import Dataset
from   py_leo_vlidort.copyvar  import _copyVar


META = ['longitude','latitude','isotime']
ASM  = ['PHIS']
MET  = ['T','PBLH','QV']
AER  = ['AIRDENS','RH','PS','DELP']
EXT  = ['ext','backscat','depol'] 
VLIDORT = ['surf_reflectance','I','ROT','rayleigh_depol_ratio','solar_zenith']
CHM  = ['O3']

AERNAMES = ['OC','SS','SU','BC','DU']
#------------------------------------ M A I N ------------------------------------

class MERGE(object):
    """
    Merge object
    """
    def __init__(self,vlidortFile,extFile,aerFile,metFile,asmFile,chmFile,outFile):
        self.vlidortFile = vlidortFile
        self.extFile = extFile
        self.aerFile = aerFile
        self.metFile = metFile
        self.asmFile = asmFile
        self.chmFile = chmFile
        self.outFile = outFile
        self.verbose = True


    def writenc(self,zlib=True):
        """
        Write a NetCDF file of merged output
        """
        km = 72

        if not os.path.exists(os.path.dirname(self.outFile)):
            os.makedirs(os.path.dirname(self.outFile))

        # Open NC file
        # ------------
        nc = Dataset(self.outFile,'w',format='NETCDF4_CLASSIC')

        # Set global attributes
        # ---------------------
        nc.title = 'Merged Inputs for Goddard lidar signal simulator'
        nc.institution = 'NASA/Goddard Space Flight Center'
        nc.source = 'Global Model and Assimilation Office'
        nc.history = ''
        nc.references = 'n/a'
        nc.contact = 'Patricia Castellanos <patricia.castellanos@nasa.gov>'
        nc.Conventions = 'CF'
        nc.vlidortFile = self.vlidortFile
        nc.extFile     = self.extFile
        nc.aerFile = self.aerFile
        nc.metFile = self.metFile
        nc.asmFile = self.asmFile
        nc.chmFile = self.chmFile
        nc.vlidort = str(VLIDORT)
        nc.ext = str(EXT)
        nc.asm = str(ASM)
        nc.met = str(MET)
        nc.aer = str(AER)
        nc.chm = str(CHM)
        nc.comment = 'VLIDORT TOA reflectance simulations limited to SZA < 80'        
     

        # Open extFile for reading
        nctrj = Dataset(self.extFile)
        ntime = len(nctrj.dimensions['time'])
     
        # Create dimensions
        # -----------------
        nt = nc.createDimension('time',ntime)
        ls = nc.createDimension('ls',19)
        nz = nc.createDimension('lev',km)
        x  = nc.createDimension('x',1)
        y  = nc.createDimension('y',1)
        ch  = nc.createDimension('channel',1)

        _copyVar(nctrj,nc,u'longitude',dtype='f4',zlib=False,verbose=self.verbose)
        _copyVar(nctrj,nc,u'latitude',dtype='f4',zlib=False,verbose=self.verbose)
        _copyVar(nctrj,nc,u'time', dtype='i4',zlib=False,verbose=self.verbose)
        _copyVar(nctrj,nc,u'lev', dtype='S1',zlib=False,verbose=self.verbose)
        _copyVar(nctrj,nc,u'isotime', dtype='S1',zlib=False,verbose=self.verbose)
        _copyVar(nctrj,nc,u'x',dtype='f4',zlib=False,verbose=self.verbose)
        _copyVar(nctrj,nc,u'y',dtype='f4',zlib=False,verbose=self.verbose)   

        # Ext variables
        for var in EXT:
            _copyVar(nctrj,nc,var,dtype='f4',zlib=False,verbose=self.verbose)

        nctrj.close()

        # Speciated
        for spc in AERNAMES:
            spcFile = self.extFile[:-3]
            spcFile = spcFile + spc + '.nc4'

            nctrj = Dataset(spcFile)
            for var in EXT:
                _copyVar(nctrj,nc,var,dtype='f4',zlib=False,verbose=self.verbose,rename=var+'_'+spc)

        # ASM
        nctrj = Dataset(self.asmFile)
        for var in ASM:
            _copyVar(nctrj,nc,var,dtype='f4',zlib=False,verbose=self.verbose)

        nctrj.close()

        #MET
        nctrj = Dataset(self.metFile)
        for var in MET:
            _copyVar(nctrj,nc,var,dtype='f4',zlib=False,verbose=self.verbose)

        nctrj.close()

        #AER
        nctrj = Dataset(self.aerFile)
        for var in AER:
            _copyVar(nctrj,nc,var,dtype='f4',zlib=False,verbose=self.verbose)

        nctrj.close()

        #CHM
        nctrj = Dataset(self.chmFile)
        for var in CHM:
            _copyVar(nctrj,nc,var,dtype='f4',zlib=False,verbose=self.verbose)

        nctrj.close() 

        #VLIDORT
        nctrj = Dataset(self.vlidortFile)
        for var in VLIDORT:
            _copyVar(nctrj,nc,var,dtype='f4',zlib=False,verbose=self.verbose)       

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
    vlidortTemplate    = cf('outDir')    + '/' + cf('outFile')
    extTemplate    = cf('outDir')    + '/' + '%instname-g5nr.lc.ext.%nymd_%hour00z.%chnm.nc4'
    outTemplate    = cf('outDir')    + '/' + '%instname-g5nr.lc2.%nymd_%hour00z.%chnm.nc4'
    channel        = int(cf('channel'))
    instname       = cf('instname')
    INSTNAME       = instname.upper()

    # Loop through dates, merging datafile
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

        inFile      = inTemplate.replace('%year',year).replace('%month',month).replace('%day',day).replace('%nymd',nymd).replace('%hour',hour).replace('%instname',instname).replace('%INSTNAME',INSTNAME)
        vlidortFile = vlidortTemplate.replace('%year',year).replace('%month',month).replace('%day',day).replace('%nymd',nymd).replace('%hour',hour).replace('%chd',get_chd(channel)).replace('%instname',instname).replace('%INSTNAME',INSTNAME)
        extFile     = extTemplate.replace('%year',year).replace('%month',month).replace('%day',day).replace('%nymd',nymd).replace('%hour',hour).replace('%ch',str(channel)).replace('%instname',instname).replace('%INSTNAME',INSTNAME)
        asmFile     = inFile.replace('%col','asm_Nx')
        metFile     = inFile.replace('%col','met_Nv')
        aerFile     = inFile.replace('%col','aer_Nv')
        chmFile     = inFile.replace('%col','chm_Nv')

        outFile     = outTemplate.replace('%year',year).replace('%month',month).replace('%day',day).replace('%nymd',nymd).replace('%hour',hour).replace('%ch',str(channel)).replace('%instname',instname).replace('%INSTNAME',INSTNAME).replace('LevelC','LevelC2')


        # Initialize MERGE class and write new outfile
        # -----------------------------------------------------------
        print '++++Merging with the following arguments+++'
        print '>>>extFile:    ',extFile
        print '>>>vlidortFile:   ',vlidortFile
        print '>>>asmFile:    ',asmFile
        print '>>>metFile:  ',metFile
        print '>>>naerile:  ',aerFile
        print '>>>chmFile:    ',chmFile
        print '>>>outFile:    ',outFile
        print '++++End of arguments+++'
        if not args.dryrun:
            merge = MERGE(vlidortFile,extFile,aerFile,metFile,asmFile,chmFile,outFile)
            merge.writenc()

        date += Dt
