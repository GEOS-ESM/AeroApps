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

from datetime        import datetime, timedelta
from dateutil.parser import parse         as isoparser

from   MAPL  import config
from   leo_vlidort import LEO_VLIDORT, get_chd, SDS_AER, SDS_CLD, SDS_INV, nMom
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


class PACE_VLIDORT(LEO_VLIDORT,LEVELBCS,PACE):
    """
    Everything needed for calling VLIDORT
    Inherit methods from LEO_VLIDORT
    GEOS-5 has already been sampled on lidar track
    """
    def __init__(self,LGinFile,LBinFile,LCinFile,outFile,rcFile,albedoType,
                channel,
                brdfFile=None,
                ndviFile=None,
                lcFile=None,
                lerFile=None,
                verbose=False,
                extOnly=False,
                distOnly=False,
                outFileDist=None,
                pool=None):

        self.SDS_AER      = SDS_AER
        self.SDS_CLD      = SDS_CLD
        self.SDS_INV      = SDS_INV
        self.SDS_GEOM     = SDS_GEOM
        self.LGinFile     = LGinFile
        self.LBinFile     = LBinFile
        self.LCinFile     = LCinFile
        self.outFile      = outFile
        self.outFileDist  = outFileDist
        self.albedoType   = albedoType
        self.rcFile       = rcFile
        self.channel      = channel
        self.verbose      = verbose
        self.nMom         = nMom
        self.brdfFile     = brdfFile
        self.lcFile       = lcFile
        self.ndviFile     = ndviFile
        self.lerFile      = lerFile

        if not extOnly:
            # initialize empty lists
            for sds in self.SDS_AER+self.SDS_CLD+self.SDS_INV:
                self.__dict__[sds] = []

            # Read in geometry and model data
            self.readGeom()
            self.readSampledGEOS()


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
                v = v.T

            self.__dict__[sds] = v  #(nobs,nz)

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
            v = v.T

            self.__dict__[sds] = v  #(nobs,nz)

            if self.verbose: 
                print 'Read ',sds

              


#------------------------------------ M A I N ------------------------------------

if __name__ == "__main__":

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
    albedoType   = 'MODIS_BRDF_BPDF'

    channel  = 470
    chd      = get_chd(channel)
    outDir    = '/nobackup/3/pcastell/POLAR_LIDAR/CALIPSO/LevelC/Y{}/M{}'.format(date.year,str(date.month).zfill(2))
    outFile   = '{}/calipso-g5nr.vlidort.vector.MCD43C.{}_{}z_{}nm.nc4'.format(outDir,nymd,hour,chd)
    
    rcFile   = 'Aod_EOS.rc'
    VZAname  = 'POLDER'
    orbit    = 'LEO'
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
                            distOnly=True)

   
    # Run ext_sampler
    # vlidort.runExt()

    # Run VLIDORT
    # if vlidort.nobs > 0:
    #     vlidort.runVLIDORT()

