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

        self.offview = ~self.offview
        self.ixs = 0
        self.ixe = 100
        self.offview[self.ixs:self.ixe,self.ixs:self.ixe] = False

        self.nobs = np.sum(~self.offview)
        self.km   = 72
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
        nc       = Dataset(AERFile)

        for sds in self.SDS_AER:
            var = nc.variables[sds]
            # v_ = np.zeros([var.shape[0],self.nobs])
            # for k in range(var.shape[0]):   
            #     v_[k,:] = var[k,:,:][~self.offview]   

            if len(var.shape) == 3:
                v_ = var[0,self.ixs:self.ixe,self.ixs:self.ixe].reshape(self.nobs) 
            else:
                v_ = var[0,:,self.ixs:self.ixe,self.ixs:self.ixe].reshape([self.km,self.nobs]) 
                v_ = v_.T

            self.__dict__[sds] = v_  #(nobs,nz)

            if self.verbose: 
                print 'Read ',sds


        if self.verbose: 
            print 'opening file',self.LCinFile
        nc       = Dataset(self.LCinFile)

        for sds in self.SDS_CLD:
            # var = nc.variables[sds]
            # v_ = np.zeros([var.shape[0],self.nobs])
            # for k in range(var.shape[0]):   
            #     v_[k,:] = var[k,:,:][~self.offview]   

            v_ = nc.variables[sds][:,self.ixs:self.ixe,self.ixs:self.ixe].reshape([self.km,self.nobs])            
            self.__dict__[sds] = v_.T  #(nobs,nz)

            if self.verbose: 
                print 'Read ',sds

        nc.close()          


#------------------------------------ M A I N ------------------------------------

if __name__ == "__main__":

    date     = datetime(2006,03,24,00,50)
    gdate    = datetime(2020,date.month,date.day,date.hour,date.minute)
    nymd     = str(date.date()).replace('-','')
    hhmm     = date.strftime('%H%M')
    hour     = date.strftime('%H')
    format   = 'NETCDF4_CLASSIC'

    LGinDir      = '/nobackup/3/pcastell/PACE/L1B/{}'.format(gdate.strftime('Y%Y/M%m/D%d'))
    LGinFile     = '{}/OCI{}00.L1B_PACE.nc'.format(LGinDir,gdate.strftime('%Y%j%H%M'))
    LBinDir      = '/nobackup/3/pcastell/PACE/LevelB/{}'.format(date.strftime('Y%Y/M%m/D%d'))
    LBinFile     = '{}/pace-g5nr.lb.%col.{}_{}00.nc4'.format(LBinDir,nymd,hhmm)
    LCinDir      = '/nobackup/3/pcastell/PACE/LevelC/{}'.format(date.strftime('Y%Y/M%m/D%d'))
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

