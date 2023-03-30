#!/usr/bin/env python3

"""
    Subset the geometry for one orbit
    closest to prime meridion

    Patricia Castellanos, Sep 2019
"""

import os
import sys
from glob import glob
import numpy as np
from scipy import signal
from netCDF4 import Dataset
from datetime import datetime, timedelta

class ORBITS(object):
    def __init__(self,rootDir,date):

        # get filename for the date
        dateDir = date.strftime('Y%Y/M%m/D%d')
        inDir = '{}/{}'.format(rootDir,dateDir) 
        filelist = glob('{}/ss450-g5nr.lb2.polar07.*'.format(inDir))

        # read in time series of data
        sza = []
        lon = []
        lat = []
        for f in filelist:
            nc = Dataset(f)
            v  = nc.variables['sza_ss'][:,4]
            sza.append(v)
            v = nc.variables['trjLon'][:]
            lon.append(v)
            v = nc.variables['trjLat'][:]
            lat.append(v)
            nc.close()

        lon = np.concatenate(lon)
        lat = np.concatenate(lat)
        sza = np.concatenate(sza)

        # get all the peaks in the time series of sza
        # these are the orbit start indeces
        a = signal.find_peaks(sza,width=10)
        a = a[0]

        # find the orbit that is close to prime meridion
        peak = None
        for i in range(len(a)-1):
            istart = a[i]
            iend   = a[i+1]
            llon = lon[istart:iend]
            llat = lat[istart:iend]
            ssza = sza[istart:iend]
            if any((llon<10) & (llon>-10) & (llat<10) & (llat>-10)& (ssza<90)):
                peak = i

        if peak is None:
            # try again
            for i in range(len(a)-1):
                istart = a[i]
                iend   = a[i+1]
                llon = lon[istart:iend]
                llat = lat[istart:iend]
                ssza = sza[istart:iend]
                if any((llon<17) & (llon>-17) & (llat<10) & (llat>-10)& (ssza<90)):
                    peak = i

        if peak is None:
            print("coulnd't find peak! Exiting")
            print(date.isoformat())
            sys.exit()
      
        self.peak   = peak
        self.a      = a 
        self.istart = a[peak]
        self.iend   = a[peak+1]
        self.filelist = filelist
        self.inDir    = inDir

#        self.lon = lon
#        self.lat = lat
#        self.sza = sza        
        # create single orbit file
        self.writeNC(date)

    def writeNC(self,date):
        nci = Dataset(self.filelist[0])


        outName = 'ss450-g5nr.orbit.polar07.{}.nc4'.format(date.strftime('%Y%m%d'))
        outFile = '{}/{}'.format(self.inDir,outName)
        print(outFile)
        nc = Dataset(outFile,'w')

        ntime = self.iend - self.istart

        # copy over global attributes
        for att in nci.ncattrs():
            nc.setncattr(att,nci.getncattr(att))     

        # set dimensions
        for dim in nci.dimensions:
            if dim == 'time':
                nc.createDimension(dim,ntime)
            else:
                nc.createDimension(dim,len(nci.dimensions[dim]))

        sds = list(nci.variables.keys())
        sds.remove('time')
        sds.remove('time_ss')
        nci.close()
        # set variables
        for vname in sds:
            
            vardata = []
            nci = Dataset(self.filelist[0])
            var = nci.variables[vname]
            dim = var.dimensions            
            if vname == 'isotime':
                varo = nc.createVariable(vname,'S1',dim,zlib=True)
            else:
                varo = nc.createVariable(vname,'f4',dim,zlib=True)
            for att in var.ncattrs():
                varo.setncattr(att,var.getncattr(att))
            if 'time' not in var.dimensions:
                #print vname
                varo[:] = var[:]
            else:
                ndim = len(var.dimensions)             
                nci.close()
  
                for fname in self.filelist:
                    nci = Dataset(fname)
                    vardata.append(nci.variables[vname][:])
                    nci.close()

                vardata = np.concatenate(vardata,axis=0)                
                #print vname,vardata.shape
                if ndim == 1:
                    vardata = vardata[self.istart:self.iend]
                elif ndim == 2:
                    vardata = vardata[self.istart:self.iend,:]
                elif ndim == 3:
                    vardata = vardata[self.istart:self.iend,:,:]

                varo[:] = vardata
               
                

        nc.close()

#------------------------------------ M A I N ------------------------------------

if __name__ == "__main__":
    rootDir   = '/discover/nobackup/projects/gmao/osse2/pub/c1440_NR/OBS/A-CCP/SS450/LevelB/'
    startdate = datetime(2006,0o7,0o3,00)
    enddate   = datetime(2007,0o1,0o1,00)
    DT        = timedelta(days=1)

    date = startdate
    while date < enddate:
        ORBITS(rootDir,date)
        date += DT
