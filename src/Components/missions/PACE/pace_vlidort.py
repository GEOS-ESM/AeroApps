#!/usr/bin/env python

"""
    PACE specific version of LEO_VLIDORT
    Calculates polarized TOA radiance for a PACE swath.
    Read in geometry from PACE L1B File
    Model fields have already been sampled. 
    Cloud simulation already created.

    Adapted from leo_vlidort.py
    Patricia Castellanos, Sep 2018

"""

import os
from   netCDF4 import Dataset
import numpy   as np

from   mieobs  import  getAOPvector, getAOPint, getEdgeVars
from datetime        import datetime, timedelta
from dateutil.parser import parse         as isoparser

from   MAPL  import config
from   leo_vlidort import LEO_VLIDORT, get_chd, SDS_AER, SDS_CLD, SDS_WIND, AERNAMES, MIENAMES
from   pace  import PACE, LEVELBCS
import multiprocessing
import READ_NC_

SDS_GEOM = {'longitude'      : 'geolocation_data',
            'sensor_azimuth' : 'geolocation_data',
            'sensor_zenith'  : 'geolocation_data',
            'solar_azimuth'  : 'geolocation_data',
            'solar_zenith'   : 'geolocation_data'}

WaterAlbedos = 'CX',

class HOLDER(PACE):
    def __init__(self):
        pass


class MieVARS(object):
    """
    container for mie vars calculations
    """
    pass


def unwrap_doMie(arg, **kwarg):
    return doMie(*arg, **kwarg)

def doMie(AERDATA,channel,rcFile,nMom,t):
    NAMES = AERNAMES + ['PS','DELP','RH','AIRDENS']
    inVars = MieVARS()
    for i,sds in enumerate(NAMES):
        var = AERDATA[i]
        if len(var.shape) ==2:
            v = var[:,t]
            v.shape = v.shape + (1,)
            inVars.__dict__[sds] = v
        elif len(var.shape) ==1:
            v = np.array(var[t])
            v.shape = (1,)
            inVars.__dict__[sds] = v

    tau,ssa,g,pmom = getAOPvector(inVars,channel,
                             vnames=AERNAMES,
                             Verbose=False,
                             rcfile=rcFile,
                             nMom=nMom)

    vol, area, refr, refi, reff = getAOPint(inVars,channel,
                             vnames=AERNAMES,
                             Verbose=False,
                             rcfile=rcFile)

    return tau,ssa,g,pmom,refr,refi,vol

class PACE_VLIDORT(LEO_VLIDORT,LEVELBCS,PACE):
    """
    Everything needed for calling VLIDORT
    Inherit methods from LEO_VLIDORT
    GEOS-5 has already been sampled on LEO track
    """
    # def __init__(self,LGinFile,LBinFile,LCinFile,outFile,rcFile,albedoType,
    #             channel,
    #             brdfFile=None,
    #             ndviFile=None,
    #             lcFile=None,
    #             lerFile=None,
    #             verbose=False,
    #             extOnly=False,
    #             distOnly=False,
    #             outFileDist=None,
    #             pool=None):

    #     self.SDS_AER      = SDS_AER
    #     self.SDS_CLD      = SDS_CLD
    #     self.SDS_INV      = SDS_INV
    #     self.SDS_GEOM     = SDS_GEOM
    #     self.LGinFile     = LGinFile
    #     self.LBinFile     = LBinFile
    #     self.LCinFile     = LCinFile
    #     self.outFile      = outFile
    #     self.outFileDist  = outFileDist
    #     self.albedoType   = albedoType
    #     self.rcFile       = rcFile
    #     self.channel      = channel
    #     self.verbose      = verbose
    #     self.nMom         = nMom
    #     self.brdfFile     = brdfFile
    #     self.lcFile       = lcFile
    #     self.ndviFile     = ndviFile
    #     self.lerFile      = lerFile

    #     if not extOnly:
    #         # initialize empty lists
    #         for sds in self.SDS_AER+self.SDS_CLD+self.SDS_INV:
    #             self.__dict__[sds] = []

    #         # Read in geometry and model data
    #         self.readGeom()
    #         self.readSampledGEOS()


    #---
    def readGeom(self):
        """
        Read in PACE geometry
        """
        if self.verbose: 
            print 'opening file',self.LGinFile

        self.geom = HOLDER()
        PACE.__init__(self.geom,self.LGinFile,SDS=SDS_GEOM)

        self.offview = self.geom.lon[0].mask
        self.nobs = np.sum(~self.offview)
        self.km   = 72
        self.offview.shape = (1,) + self.offview.shape 
        self.offview3d = np.repeat(self.offview,self.km,axis=0)
        self.offview.shape = self.offview.shape[1:]
        self.shape = self.geom.lon[0].shape

        # Store flattened arrays
        for sds in self.SDS_GEOM:
            self.__dict__[sds] =  self.geom.__dict__[sds][0][~self.offview]

        self.geom = None

    #---
    def readSampledGEOS(self):
        """
        Read in model sampled track
        """
        col = 'aer_Nv'
        AERFile = self.LBinFile.replace('%col',col)
        if self.verbose: 
            print 'opening file',AERFile
 
        for sds in self.SDS_AER:
            if sds == 'PS':
                var = np.zeros([self.shape[1],self.shape[0],1],dtype=np.float32,order='F')
                READ_NC_.py_readvar3d(sds,AERFile,var)
                v = np.squeeze(var).T
                v = v[~self.offview]

            else:
                var = np.zeros([self.shape[1],self.shape[0],self.km,1],dtype=np.float32,order='F')
                READ_NC_.py_readvar4d(sds,AERFile,var)

                v = np.squeeze(var).T
                v = v[~self.offview3d]
                v = v.reshape([self.km,self.nobs])

            self.__dict__[sds] = v  #(nz,nobs)

            if self.verbose: 
                print 'Read ',sds


        if self.verbose: 
            print 'opening file',self.LCinFile

        for sds in self.SDS_CLD:
            var = np.zeros([self.shape[1],self.shape[0],self.km],dtype=np.float32,order='F')
            READ_NC_.py_readvar3d(sds,LCinFile,var)

            v = np.squeeze(var).T
            v = v[~self.offview3d]
            v = v.reshape([self.km,self.nobs])

            self.__dict__[sds] = v  #(nz,nobs)

            if self.verbose: 
                print 'Read ',sds


    #---
    def read10mWind(self):
        """
        Read in model sampled track wind data
        """
        col = 'met_Nv'
        self.SDS_WIND = SDS_WIND
        METFile = self.LBinFile.replace('%col',col)
        if self.verbose: 
            print 'opening file',METFile
 
        for sds in self.SDS_WIND:
            var = np.zeros([self.shape[1],self.shape[0],1],dtype=np.float32,order='F')
            READ_NC_.py_readvar3d(sds,METFile,var)
            v = np.squeeze(var).T
            v = v[~self.offview]

            self.__dict__[sds] = v  #(nobs,nz)

            if self.verbose: 
                print 'Read ',sds
              


#------------------------------------ M A I N ------------------------------------

if __name__ == "__main__":
    # # start fork here so large data is not copied to children
    nproc = int(multiprocessing.cpu_count()*0.5)    
    pool = multiprocessing.Pool(nproc) 

    date     = datetime(2006,03,24,00,50)
    gdate    = datetime(2020,date.month,date.day,date.hour,date.minute)
    nymd     = str(date.date()).replace('-','')
    hhmm     = date.strftime('%H%M')
    hour     = date.strftime('%H')
    format   = 'NETCDF4_CLASSIC'

    LGinDir      = '/nobackup/PACE/L1B/{}'.format(gdate.strftime('Y%Y/M%m/D%d'))
    LGinFile     = '{}/OCI{}00.L1B_PACE.nc'.format(LGinDir,gdate.strftime('%Y%j%H%M'))
    LBinDir      = '/nobackup/PACE/LevelB/{}'.format(date.strftime('Y%Y/M%m/D%d'))
    LBinFile     = '{}/pace-g5nr.lb.%col.{}_{}00.nc4'.format(LBinDir,nymd,hhmm)
    LCinDir      = '/nobackup/PACE/LevelC/{}'.format(date.strftime('Y%Y/M%m/D%d'))
    LCinFile     = '{}/pace-g5nr.TOTWPDF-GCOP-SKEWT.{}_{}00.nc4'.format(LCinDir,nymd,hhmm)


    brdfDir      = '/nobackup/3/pcastell/POLAR_LIDAR/CALIPSO/BRDF/MCD43C1/006/Y{}/M{}'.format(date.year,str(date.month).zfill(2))
    brdfFile     = '{}/calipso-g5nr.lb2.brdf.{}_{}z.nc4'.format(brdfDir,nymd,hour)
    ndviDir      = '/nobackup/3/pcastell/POLAR_LIDAR/CALIPSO/BPDF/NDVI/MYD13C2/006/Y{}/M{}'.format(date.year,str(date.month).zfill(2))
    ndviFile     = '{}/calipso-g5nr.lb2.ndvi.{}_{}z.nc4'.format(ndviDir,nymd,hour)
    lcDir        = '/nobackup/3/pcastell/POLAR_LIDAR/CALIPSO/BPDF/LAND_COVER/MCD12C1/051/Y{}/M{}'.format(date.year,str(date.month).zfill(2))
    lcFile       = '{}/calipso-g5nr.lb2.land_cover.{}_{}z.nc4'.format(lcDir,nymd,hour)    
    albedoType   = 'CX'

    channel  = 470
    chd      = get_chd(channel)
    outDir    = '/nobackup/PACE/LevelC/Y{}/M{}'.format(date.year,str(date.month).zfill(2))
    outFile   = '{}/pace-g5nr.vlidort.vector.CX.{}_{}_{}nm.nc4'.format(outDir,nymd,hhmm,chd)
    
    rcFile   = 'Aod_EOS.rc'
    verbose  = True

    # Initialize VLIDORT class getting aerosol optical properties
    # -----------------------------------------------------------
    v  = PACE_VLIDORT(LGinFile,LBinFile,LCinFile,outFile,rcFile,
                            albedoType,
                            channel,
                            brdfFile=brdfFile,
                            ndviFile=ndviFile,
                            lcFile=lcFile,
                            verbose=verbose,
                            distOnly=False)

    if not v.extOnly:
        # Calculate aerosol optical properties

        #doMie(AERDATA,channel,rcFile,nMom,t)

        NAMES = AERNAMES + ['PS','DELP','RH','AIRDENS']
        nbatch = 10


        #for sii in np.arange(0,v.nobs,nbatch): 
        for sii in np.arange(1):   
            eii = sii + nbatch 
            nobs = nbatch
            if eii > v.nobs:
                eii = v.nobs
                nobs = eii - sii

            # for sds in MIENAMES:
            #     if sds == 'pmom':
            #         v.__dict__[sds] = np.zeros([v.km,1,nobs,v.nMom,4])
            #     else:
            #         v.__dict__[sds] = np.zeros([v.km,1,nobs]) #(km,nch,nobs)

            
            args = HOLDER()
            argsList = []
            for sds in NAMES:
                var = v.__dict__[sds]
            
                if len(var.shape) == 2:
                    args.__dict__[sds] = var[:,sii:eii]
                else:
                    args.__dict__[sds] = var[sii:eii]

                argsList.append(args.__dict__[sds])

            argsALL = zip(nobs*[argsList],nobs*[v.channel],nobs*[v.rcFile],nobs*[v.nMom],range(nobs))
            result = pool.map(unwrap_doMie,argsALL)
            # for r in result:
            #     # tau, ssa, g, pmom,refi,refr,vol = r
            #     for sds,miedata in zip(MIENAMES,r):
            #         v.__dict__[sds] = miedata
            #         # if sds == 'pmom':
            #         #     v.__dict__[sds][:,:,sii:eii,:,:] = miedata
            #         # else:
            #         #     v.__dict__[sds][:,:,sii:eii] = miedata

   
    # Run ext_sampler
    # vlidort.runExt()

    # Run VLIDORT
    # if vlidort.nobs > 0:
    #     vlidort.runVLIDORT()

    pool.close()
    pool.join()

