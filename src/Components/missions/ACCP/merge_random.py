#!/usr/bin/env python3

"""
    Puts all needed lidar signal simulator variables in one file

    Patricia Castellanos, Dec, 2019

"""

import os
import argparse
from   datetime        import datetime, timedelta
from   dateutil.parser import parse         as isoparser
from   dateutil.relativedelta import relativedelta
from   MAPL.config      import Config
from   py_leo_vlidort.lidar_vlidort   import get_chd
import numpy  as np
from   netCDF4 import Dataset
from   py_leo_vlidort.copyvar  import _copyVar


AERNAMES = ['OC','SS','SU','BC','DU']
#------------------------------------ M A I N ------------------------------------

class MERGE(object):
    """
    Merge object
    """
    def __init__(self,date,randomFile,args):
        self.randomFile = randomFile
        self.args       = args

        # date
        nymd  = str(date.date()).replace('-','')
        year  = str(date.year)
        month = str(date.month).zfill(2)


        # get random draws
        self.read_random()

        # --- orbitname ---
        cf              = Config(args.orbit_pcf,delim=' = ')
        orbitname       = cf('orbitname')
        ORBITNAME       = orbitname.upper()

        # --- instrument name ---
        cf              = Config(args.inst_pcf,delim=' = ')
        instname        = cf('instname')

        # --- instrument name ---
        cf              = Config(args.ext_pcf,delim=' = ')
        channels        = cf('CHANNELS')
        channels        = channels.replace(' ','')
        if ',' in channels:
            channels = channels.split(',')
        else:
            channels = [channels]
        species         = AERNAMES

        # -- inFiles (Level B) ---
        cf             = Config(args.track_pcf,delim=' = ')
        inTemplate     = cf('inDir')     + '/' + cf('inFile')
        outDir         = cf('inDir').split('/')[:-1]
        outDir         = '/'.join(outDir)
        outTemplate    = outDir + '/' + cf('inFile')

        # loop through collections
        cols    = 'aer_Nv','chm_Nv','asm_Nx','met_Nv',instname,'aer_SD'
#        cols    = []
        inTemplate  = inTemplate.replace('%orbitname',orbitname).replace('%ORBITNAME',ORBITNAME)
        outTemplate = outTemplate.replace('%orbitname',orbitname).replace('%ORBITNAME',ORBITNAME)
        for col in cols:
            print('Working on ',col)
            inFile = inTemplate.replace('%col',col)
            outFile = outTemplate.replace('%col',col+'.random')

            # create outfile
            day = '01'
            hour = '01'
            inFiled = inFile.replace('%year',year).replace('%month',month).replace('%day',day).replace('%nymd',nymd).replace('%hour',hour)
            hour = '00'
            outFile = outFile.replace('%year',year).replace('%month',month).replace('%nymd',nymd).replace('%hour',hour)
            
            self.make_outFile(inFiled,outFile)
            self.copyRandom(inFile,outFile) 


        # --- inFiles (Level B Surface) ---
        cols = 'brdf','ndvi','lc','ler'
#        cols = []
        for col in cols:
            try:
                inTemplate = cf(col+'Dir') + '/' + cf(col+'File')
                outDir       = cf(col+'Dir').split('/')[:-1]
                outDir         = '/'.join(outDir)
                outTemplate    = outDir + '/' + cf(col+'File')
            except:
                inTemplate = None
        
            if inTemplate is not None:
                print('Working on ',col)
                inFile  = inTemplate.replace('%orbitname',orbitname).replace('%ORBITNAME',ORBITNAME)
                outFile = outTemplate.replace('%orbitname',orbitname).replace('%ORBITNAME',ORBITNAME)
            
                outFile = outFile.replace('.'+col+'.','.'+col+'.random.')
                # create outfile
                day = '01'
                hour = '01'
                inFiled = inFile.replace('%year',year).replace('%month',month).replace('%day',day).replace('%nymd',nymd).replace('%hour',hour)
                hour = '00'
                outFile = outFile.replace('%year',year).replace('%month',month).replace('%nymd',nymd).replace('%hour',hour)

                self.make_outFile(inFiled,outFile)
                self.copyRandom(inFile,outFile)


        # -- extinction files (LevelC) ---
        # total 
        inTemplate    = cf('outDir')    + '/' + '%orbitname-g5nr.lc.ext.%nymd_%hour00z.%chnm.nc4'
        outDir       = cf('outDir').split('/')[:-1]
        outDir         = '/'.join(outDir)
        outTemplate    = outDir + '/' + '%orbitname-g5nr.lc.ext.random.%nymd_%hour00z.%chnm.nc4'

        inTemplate  = inTemplate.replace('%orbitname',orbitname).replace('%ORBITNAME',ORBITNAME)
        outTemplate = outTemplate.replace('%orbitname',orbitname).replace('%ORBITNAME',ORBITNAME)
        for ch in channels:
            print('Working on ext ',ch)
            inFile = inTemplate.replace('%ch',ch)
            outFile = outTemplate.replace('%ch',ch)

            # create outfile
            day = '01'
            hour = '01'
            inFiled = inFile.replace('%year',year).replace('%month',month).replace('%day',day).replace('%nymd',nymd).replace('%hour',hour)
            hour = '00'
            outFile = outFile.replace('%year',year).replace('%month',month).replace('%nymd',nymd).replace('%hour',hour)

            self.make_outFile(inFiled,outFile)
            self.copyRandom(inFile,outFile)

        # total speciated
        inTemplate    = cf('outDir')    + '/' + '%orbitname-g5nr.lc.ext.%nymd_%hour00z.%chnm.%spc.nc4'
        outDir       = cf('outDir').split('/')[:-1]
        outDir         = '/'.join(outDir)
        outTemplate    = outDir + '/' + '%orbitname-g5nr.lc.ext.random.%nymd_%hour00z.%chnm.%spc.nc4'

        inTemplate  = inTemplate.replace('%orbitname',orbitname).replace('%ORBITNAME',ORBITNAME)
        outTemplate = outTemplate.replace('%orbitname',orbitname).replace('%ORBITNAME',ORBITNAME)
        for spc in species:
            for ch in channels:
                print('Working on ext ',ch,' ',spc)
                inFile = inTemplate.replace('%ch',ch).replace('%spc',spc)
                outFile = outTemplate.replace('%ch',ch).replace('%spc',spc)

                # create outfile
                day = '01'
                hour = '01'
                inFiled = inFile.replace('%year',year).replace('%month',month).replace('%day',day).replace('%nymd',nymd).replace('%hour',hour)
                hour = '00'
                outFile = outFile.replace('%year',year).replace('%month',month).replace('%nymd',nymd).replace('%hour',hour)

                self.make_outFile(inFiled,outFile)
                self.copyRandom(inFile,outFile)


        # fine mode
        inTemplate    = cf('outDir')    + '/' + '%orbitname-g5nr.lc.ext.finemode.%nymd_%hour00z.%chnm.nc4'
        outDir       = cf('outDir').split('/')[:-1]
        outDir         = '/'.join(outDir)
        outTemplate    = outDir + '/' + '%orbitname-g5nr.lc.ext.finemode.random.%nymd_%hour00z.%chnm.nc4'

        inTemplate  = inTemplate.replace('%orbitname',orbitname).replace('%ORBITNAME',ORBITNAME)
        outTemplate = outTemplate.replace('%orbitname',orbitname).replace('%ORBITNAME',ORBITNAME)
        for ch in channels:
            print('Working on extfine ',ch)
            inFile = inTemplate.replace('%ch',ch)
            outFile = outTemplate.replace('%ch',ch)

            # create outfile
            day = '01'
            hour = '01'
            inFiled = inFile.replace('%year',year).replace('%month',month).replace('%day',day).replace('%nymd',nymd).replace('%hour',hour)
            hour = '00'
            outFile = outFile.replace('%year',year).replace('%month',month).replace('%nymd',nymd).replace('%hour',hour)

            self.make_outFile(inFiled,outFile)
            self.copyRandom(inFile,outFile)


        # finemode speciated
        inTemplate    = cf('outDir')    + '/' + '%orbitname-g5nr.lc.ext.finemode.%nymd_%hour00z.%chnm.%spc.nc4'
        outDir       = cf('outDir').split('/')[:-1]
        outDir         = '/'.join(outDir)
        outTemplate    = outDir + '/' + '%orbitname-g5nr.lc.ext.finemode.random.%nymd_%hour00z.%chnm.%spc.nc4'

        inTemplate  = inTemplate.replace('%orbitname',orbitname).replace('%ORBITNAME',ORBITNAME)
        outTemplate = outTemplate.replace('%orbitname',orbitname).replace('%ORBITNAME',ORBITNAME)
        for spc in species:
            for ch in channels:
                print('Working on extfine ',ch,' ',spc)
                inFile = inTemplate.replace('%ch',ch).replace('%spc',spc)
                outFile = outTemplate.replace('%ch',ch).replace('%spc',spc)

                # create outfile
                day = '01'
                hour = '01'
                inFiled = inFile.replace('%year',year).replace('%month',month).replace('%day',day).replace('%nymd',nymd).replace('%hour',hour)
                hour = '00'
                outFile = outFile.replace('%year',year).replace('%month',month).replace('%nymd',nymd).replace('%hour',hour)

                self.make_outFile(inFiled,outFile)
                self.copyRandom(inFile,outFile)


        # lidar merged files
        inTemplate  = cf('outDir')    + '/' + '%orbitname-lidar-g5nr.lc2.%nymd_%hour00z.%chnm.nc4'
        outDir       = cf('outDir').split('/')[:-1]
        outDir         = '/'.join(outDir)
        outTemplate    = outDir + '/' + '%orbitname-lidar-g5nr.lc2.random.%nymd_%hour00z.%chnm.nc4'
        inTemplate = inTemplate.replace('LevelC','LevelC2')
        outTemplate = outTemplate.replace('LevelC','LevelC2')

        inTemplate  = inTemplate.replace('%orbitname',orbitname).replace('%ORBITNAME',ORBITNAME)
        outTemplate = outTemplate.replace('%orbitname',orbitname).replace('%ORBITNAME',ORBITNAME)
        channels = '355','532','1064'
        for ch in channels:
            print('working on merged ',ch)
            inFile = inTemplate.replace('%ch',ch)
            outFile = outTemplate.replace('%ch',ch)

            # create outfile
            day = '01'
            hour = '01'
            inFiled = inFile.replace('%year',year).replace('%month',month).replace('%day',day).replace('%nymd',nymd).replace('%hour',hour)
            hour = '00'
            outFile = outFile.replace('%year',year).replace('%month',month).replace('%nymd',nymd).replace('%hour',hour)

            self.make_outFile(inFiled,outFile)
            self.copyRandom(inFile,outFile)


    def copyRandom(self,inFile,outFile):
        """
        Get random pixel and write to outfile
        """

        # outFile
        nco = Dataset(outFile,'r+')

        # loop through variables with a time dimension
        for varname in nco.variables:
            print('var ',varname)
            var = nco.variables[varname]
            dim = var.dimensions
            if 'time' in dim:
                varout = []
                idim = np.argmax(np.array(dim) == 'time')
                ndim = len(dim)
                # loop through random files
                for irandom,dd in enumerate(self.random_times):
                    index = self.random_secs[irandom]
                    Rnymd,Rhh = dd.split('_')
                    Ryear     = Rnymd[0:4]
                    Rmonth    = Rnymd[4:6]
                    Rday      = Rnymd[6:]
                    Rhour     = Rhh[0:2]
                    inFiled = inFile.replace('%year',Ryear).replace('%month',Rmonth).replace('%day',Rday).replace('%nymd',Rnymd).replace('%hour',Rhour)


                    #inFile
                    nci = Dataset(inFiled)

                    if ndim == 1:
                        data = nci.variables[varname][index:index+1]
                    else:
                        if ndim == 2:
                            if idim == 0:
                                data = nci.variables[varname][index:index+1,:]
                            else:
                                data = nci.variables[varname][:,index:index+1]
                        elif ndim == 3:
                            if idim == 0:
                                data = nci.variables[varname][index:index+1,:,:]
                            elif idim == 1:
                                data = nci.variables[varname][:,index:index+1,:]
                            else:
                                data = nci.variables[varname][:,:,index:index+1]
                        elif ndim == 4:
                            if idim == 0:
                                data = nci.variables[varname][index:index+1,:,:,:]
                            elif idim == 1:
                                data = nci.variables[varname][:,index:index+1,:,:]
                            elif idim == 2:
                                data = nci.variables[varname][:,:,index:index+1,:]
                            else:
                                data = nci.variables[varname][:,:,:,index:index+1]
                        else:
                            print('Dimension error.  Add another dimension')
                    varout.append(data)
                    nci.close()
             
                varout = np.concatenate(varout,axis=idim)                 
                var[:] = varout
        nco.close()    
    def make_outFile(self,inFile,outFile):
        """
        Make an outfile with same dimensions and attributes as inFile
        Time dimension is unlimited so it can be appended to
        """

        # inFile
        nci = Dataset(inFile)

        # outFile
        nco = Dataset(outFile,'w',format='NETCDF4_CLASSIC')

        # file attributes
        for atr in nci.__dict__:
            text = nci.__dict__[atr]
            nco.setncattr(atr,text)

        # dimensions
        for dim in nci.dimensions:
            if dim == 'time':
                nco.createDimension(dim,self.nrandom)
            else:
                nco.createDimension(dim,len(nci.dimensions[dim])) 
        
        # variables
        for varname in nci.variables:
            var = nci.variables[varname]
            dim = var.dimensions
            dtype = var.dtype
            obj = nco.createVariable(varname,dtype,dimensions=dim,zlib=True)
            for atr in var.__dict__:
                text = var.__dict__[atr]
                obj.setncattr(atr,text)

            # if no time dimension, write it out
            if 'time' not in dim:
                data = var[:]
                obj[:] = data
            
        nci.close()
        nco.close()


    def read_random(self):
        """
        Read dates and times from random draw file
        """
        f = open(self.randomFile)
        ndraw, nang = f.readline().split()
        ndraw = int(ndraw)
        nang  = int(nang)

        lines = []
        for ll in f:
            lines.append(ll)
        f.close()

        times = []
        secs  = []
        for ll in np.arange(0,ndraw*nang,nang):
            lon,lat,time,sec,dsec,ang1,ang2,ang3,ang4,ang5 = lines[ll].split(',')
            times.append(time)
            secs.append(sec)

        self.random_times = times
        self.random_secs  = np.array(secs).astype(int)
        self.nrandom      = ndraw


if __name__ == "__main__":

#   Defaults
#   ---------
    DT_mons = 1
   
#   Parse command line options
#   --------------------------

    parser = argparse.ArgumentParser()
    parser.add_argument("iso_t1",
                        help="starting iso time")
    parser.add_argument("iso_t2",
                        help="ending iso time")

    parser.add_argument("track_pcf",
                        help="prep config file with with track input file names")

    parser.add_argument("orbit_pcf",
                        help="prep config file with orbit variables")

    parser.add_argument("inst_pcf",
                        help="prep config file with instrument variables")

    parser.add_argument("ext_pcf",
                        help="prep config file with extinction sampling variables")

    parser.add_argument("-D","--DT_mons", default=DT_mons, type=int,
                        help="Timestep in months for each file (default=%i)"%DT_mons)

    parser.add_argument("-v", "--verbose",action="store_true",
                        help="Verbose mode (default=False).")

    parser.add_argument("-r", "--dryrun",action="store_true",
                        help="do a dry run (default=False).")    

    args = parser.parse_args()


    # Parse prep config
    # -----------------
    
    # --- random draws ---
    cf             = Config(args.track_pcf,delim=' = ')
    LbDir          =  cf('inDir')
    rootDir        = LbDir.split('/')[:-4]
    rootDir        = '/'.join(rootDir)
    randomTemplate = rootDir + '/random_draw/%orbitname.random.%year%month.txt'

    # --- orbitname ---
    cf              = Config(args.orbit_pcf,delim=' = ')
    orbitname       = cf('orbitname')
    ORBITNAME       = orbitname.upper()



    # Loop through dates, merging datafile
    # ------------------------------------
    date      = isoparser(args.iso_t1)
    enddate   = isoparser(args.iso_t2)
    Dt        = relativedelta(months=args.DT_mons)

    while date < enddate:
        year  = str(date.year)
        month = str(date.month).zfill(2)

        randomFile = randomTemplate.replace('%ORBITNAME',ORBITNAME).replace('%orbitname',orbitname).replace('%year',year).replace('%month',month)


        # Initialize MERGE class and write new outfile
        # -----------------------------------------------------------
        print('++++Merging from random draws+++')
        print('>>>randomFile:    ',randomFile)
        print('++++End of arguments+++')
        if not args.dryrun:
            merge = MERGE(date,randomFile,args)
#            merge.writenc()


        date += Dt
