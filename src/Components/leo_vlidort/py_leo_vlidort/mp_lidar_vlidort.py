#!/usr/bin/env python

"""
    Calculates polarized TOA radiance for a lidar track.
    Model fields have already been sampled using trj_sampler

    Adapted from ext_sampler.py
    Adapted from polar_vlidort.py
    Patricia Castellanos, Dec, 2019

"""

import os
from   netCDF4 import Dataset
from   mieobs  import  getAOPvector, getEdgeVars
import numpy   as np
from py_leo_vlidort.vlidort import VLIDORT, get_chd, WrapperFuncs, CX_run, LAMBERTIAN_run, MODIS_BRDF_run
from datetime        import datetime, timedelta
from dateutil.parser import parse         as isoparser

from MAPL.constants import *
from py_leo_vlidort import LidarAngles_    
 
#from copyvar  import _copyVar
from scipy import interpolate
from   MAPL  import config
from multiprocessing import Pool

# Generic Lists of Varnames and Units
VNAMES_DU = ['DU001','DU002','DU003','DU004','DU005']
VNAMES_SS = ['SS001','SS002','SS003','SS004','SS005']
VNAMES_BC = ['BCPHOBIC','BCPHILIC']
VNAMES_OC = ['OCPHOBIC','OCPHILIC']
VNAMES_SU = ['SO4']

META    = ['DELP','PS','RH','AIRDENS','LONGITUDE','LATITUDE','isotime']
AERNAMES = VNAMES_SU + VNAMES_SS + VNAMES_OC + VNAMES_BC + VNAMES_DU
SDS_AER = META + AERNAMES
SDS_MET = [] #[CLDTOT]
SDS_INV = ['FRLAND']
SDS_CX = ['U10M','V10M']

ncALIAS = {'LONGITUDE': 'trjLon',
           'LATITUDE': 'trjLat'}

nMom     = 300

SurfaceFuncs = {'MODIS_BRDF'     : 'readSampledMODISBRDF',
                'MODIS_BRDF_BPDF': 'readSampledMODISBRDF',
                'LAMBERTIAN'     : 'readSampledLER',
                'CX'             : 'readSampledWindCX'}

MISSING = -1.e+20

class LIDAR_VLIDORT(VLIDORT):
    """
    Everything needed for calling VLIDORT
    GEOS-5 has already been sampled on lidar track
    """
    def __init__(self,inFile,outFile,rcFile,albedoType,
                channel,hgtss,
                nstreams=12,
                plane_parallel=True,
                brdfFile=None,
                ndviFile=None,
                lcFile=None,
                lerFile=None,
                verbose=False):
        self.SDS_AER     = SDS_AER
        self.SDS_MET     = SDS_MET
        self.SDS_INV     = SDS_INV
        self.SDS_CX      = SDS_CX
        self.AERNAMES    = AERNAMES
        self.inFile      = inFile
        self.outFile     = outFile
        self.albedoType  = albedoType
        self.rcFile      = rcFile
        self.channel     = channel
        self.verbose     = verbose
        self.nMom        = nMom
        self.brdfFile    = brdfFile
        self.lcFile      = lcFile
        self.ndviFile    = ndviFile
        self.lerFile     = lerFile
        self.nstreams    = nstreams
        self.plane_parallel = plane_parallel

        # initialize empty lists
        for sds in self.SDS_AER+self.SDS_MET+self.SDS_INV+self.SDS_CX:
            self.__dict__[sds] = []

        # Read in model data
        self.readSampledGEOS()

        # Make lists into arrays
        for sds in self.SDS_AER+self.SDS_MET:
            self.__dict__[sds] = np.concatenate(self.__dict__[sds])

        # convert isotime to datetime
        self.tyme = []
        for isotime in self.isotime:
            self.tyme.append(isoparser(''.join(isotime)))

        self.tyme = np.array(self.tyme)

        # Start out with all good obs
        self.nobs  = len(self.tyme)
        self.iGood = np.ones([self.nobs]).astype(bool)
        self.nobsLand  = len(self.tyme)
        self.iGoodLand = np.ones([self.nobs]).astype(bool)
        
        # Read in surface data
        # Intensity
        if (self.channel < 470) & ("MODIS_BRDF" in albedoType):
            albedoReader = getattr(self,'readHybridMODISBRDF')
        else:
            albedoReader = getattr(self,SurfaceFuncs[albedoType])
        albedoReader()

        # Polarization
        if 'BPDF' in albedoType:
            self.BPDFinputs()

        # Ocean
        #albedoReader = getattr(self,SurfaceFuncs['CX'])
        #albedoReader()

        # Calculate aerosol optical properties
        self.computeMie()

        # Calculate atmospheric profile properties needed for Rayleigh calc
        self.computeAtmos()

        # Calculate Scene Geometry
        # limit iGood to sza < 80
        self.VZA = 0.0
        self.VAA = 0.0
        self.hgtss = hgtss
        self.calcAngles()

        if self.nobs > 0:
            # Land-Sea Mask
            self.LandSeaMask()

    # --
    def calcAngles(self):
        SZA   = []
        SAA   = []

        for i,tyme in enumerate(self.tyme):
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



        self.SZA = np.array(SZA)
        self.SAA = np.array(SAA)

        # define RAA according to photon travel direction
        saa = self.SAA + 180.0
        I = saa >= 360.
        saa[I] = saa[I] - 360.

        RAA = self.VAA - saa
        RAA[RAA < 0] = RAA[RAA<0]+360.0
        self.RAA = RAA

        # Limit SZAs
        self.iGood = self.iGood & (self.SZA < 80)
        self.nobs = np.sum(self.iGood)     

    def runVLIDORT(self):
        """
        Calls VLIDORT 
        """

        # Only do good obs
        surfList = []
        if self.nobsLand > 0:
            self.iLand = np.arange(len(self.iGood))[self.iGood & self.iLand]
            surfList.append('Land')
        if self.nobsSea > 0:
            self.iSea = np.arange(len(self.iGood))[self.iGood & self.iSea]
            surfList.append('Sea')
        self.iGood = np.arange(len(self.iGood))[self.iGood]

        # Initiate output arrays
        ntime   = len(self.tyme)
        nlev    = self.tau.shape[0]
        self.I = np.ones([ntime])*MISSING
        self.Q = np.ones([ntime])*MISSING
        self.U = np.ones([ntime])*MISSING
        self.reflectance = np.ones([ntime])*MISSING
        self.surf_reflectance = np.ones([ntime])*MISSING
        self.BR_Q = np.ones([ntime])*MISSING
        self.BR_U = np.ones([ntime])*MISSING
        self.ROT = np.ones([ntime,nlev])*MISSING

        # Calculate ROT
        args = [self.channel, self.pe, self.ze, self.te, MISSING, self.verbose]
        vlidortWrapper = WrapperFuncs['ROT_CALC']
        ROT, depol_ratio, rc = vlidortWrapper(*args)  
        #nlev,ntime,nch
        self.ROT = np.squeeze(ROT).T 
        self.depol_ratio = depol_ratio   

        p = Pool(27) 
        # loop though LAND and SEA
        for surface in surfList:
            print 'Working on ', surface
            iGood = self.__dict__['i'+surface]
            nobs = len(iGood)
            sza  = self.SZA[iGood]
            raa  = self.RAA[iGood] 
            if surface == 'Land':
                vza = np.zeros(self.nobsLand)
            else:
                vza = np.zeros(self.nobsSea)

            tau = self.tau[:,:,iGood]
            ssa = self.ssa[:,:,iGood]
            pmom = self.pmom[:,:,iGood,:,:]
            pe   = self.pe[:,iGood]
            ze   = self.ze[:,iGood]
            te   = self.te[:,iGood]
            rot  = ROT[:,iGood,:]

            if surface == 'Land':        
                albedoType = self.albedoType

            else:
                albedoType = 'CX'

            # Get surface data
            if albedoType == 'MODIS_BRDF':
                param     = self.RTLSparam[:,:,iGood]
                kernel_wt = self.kernel_wt[:,:,iGood]
            elif albedoType == 'LAMBERTIAN':
                albedo = self.albedo[iGood,:]
            elif albedoType == 'CX':
                u10m = self.U10M[iGood]
                v10m = self.V10M[iGood]

            I = []
            Q = []
            U = []
            reflectance = []
            surf_reflectance = []
            BR_Q = []
            BR_U = []

            # run VLIDORT
            if albedoType == 'MODIS_BRDF':
                args = [(self.channel, self.nstreams, self.plane_parallel, ROT[:,i:i+1,:], depol_ratio,
                        tau[:,:,i:i+1], ssa[:,:,i:i+1], pmom[:,:,i:i+1,:,:],
                        pe[:,i:i+1], ze[:,i:i+1], te[:,i:i+1],
                        kernel_wt[:,:,i:i+1], param[:,:,i:i+1],
                        sza[i:i+1], raa[i:i+1], vza[i:i+1],
                        MISSING,
                        self.verbose) for i in range(nobs)]
                result = p.map(MODIS_BRDF_run,args)
                for r in result:
                    I_r,Q_r,U_r,reflectance_r,surf_reflectance_r,BR_Q_r,BR_U_r = r
                    I.append(I_r)
                    Q.append(Q_r)
                    U.append(U_r)
                    reflectance.append(reflectance_r)
                    surf_reflectance.append(surf_reflectance_r)
                    BR_Q.append(BR_Q_r)
                    BR_U.append(BR_U_r)
                I = np.concatenate(I)
                Q = np.concatenate(Q)
                U = np.concatenate(U)
                reflectance = np.concatenate(reflectance)
                surf_reflectance = np.concatenate(surf_reflectance)
                BR_Q = np.concatenate(BR_Q)
                BR_U = np.concatenate(BR_U)
            elif albedoType == 'LAMBERTIAN':
                args = [(self.channel, self.nstreams, self.plane_parallel, ROT[:,i:i+1,:], depol_ratio,
                        tau[:,:,i:i+1], ssa[:,:,i:i+1], pmom[:,:,i:i+1,:,:],
                        pe[:,i:i+1], ze[:,i:i+1], te[:,i:i+1],
                        albedo[i:i+1,:],
                        sza[i:i+1], raa[i:i+1], vza[i:i+1],
                        MISSING,
                        self.verbose) for i in range(nobs)]
                result = p.map(LAMBERTIAN_run,args)
                for r in result:
                    I_r,Q_r,U_r,reflectance_r,surf_reflectance_r,BR_Q_r,BR_U_r = r
                    I.append(I_r)
                    Q.append(Q_r)
                    U.append(U_r)
                    reflectance.append(reflectance_r)
                surf_reflectance = albedo
                BR_Q = np.ones(nobs)*MISSING
                BR_U = np.ones(nobs)*MISSING
                I = np.concatenate(I)
                Q = np.concatenate(Q)
                U = np.concatenate(U)
                reflectance = np.concatenate(reflectance)

            elif albedoType == 'CX':
                args = [(self.channel, self.nstreams, self.plane_parallel, ROT[:,i:i+1,:], depol_ratio,
                        tau[:,:,i:i+1], ssa[:,:,i:i+1], pmom[:,:,i:i+1,:,:],
                        pe[:,i:i+1], ze[:,i:i+1], te[:,i:i+1],
                        u10m[i:i+1], v10m[i:i+1], self.mr,
                        sza[i:i+1], raa[i:i+1], vza[i:i+1],
                        MISSING,
                        self.verbose) for i in range(nobs)]
                result = p.map(CX_run,args)
                for r in result:
                    I_r,Q_r,U_r,reflectance_r,surf_reflectance_r,BR_Q_r,BR_U_r = r
                    I.append(I_r)
                    Q.append(Q_r)
                    U.append(U_r)
                    reflectance.append(reflectance_r)
                    surf_reflectance.append(surf_reflectance_r)
                    BR_Q.append(BR_Q_r)
                    BR_U.append(BR_U_r)
                I = np.concatenate(I)
                Q = np.concatenate(Q)
                U = np.concatenate(U)
                reflectance = np.concatenate(reflectance)
                surf_reflectance = np.concatenate(surf_reflectance)
                BR_Q = np.concatenate(BR_Q)
                BR_U = np.concatenate(BR_U)

                                                
            self.I[iGood] = np.squeeze(I)
            self.reflectance[iGood] = np.squeeze(reflectance)
            self.surf_reflectance[iGood] = np.squeeze(surf_reflectance)
            self.Q[iGood] = np.squeeze(Q)
            self.U[iGood] = np.squeeze(U) 
            self.BR_Q[iGood] = np.squeeze(BR_Q)
            self.BR_U[iGood] = np.squeeze(BR_U) 


        self.writeNC()
    #---
    def writeNC (self,zlib=True):
        """
        Write a NetCDF file vlidort output
        """
        km = 72

        if not os.path.exists(os.path.dirname(self.outFile)):
            os.makedirs(os.path.dirname(self.outFile))

        # Open NC file
        # ------------
        nc = Dataset(self.outFile,'w',format='NETCDF4_CLASSIC')

        # Set global attributes
        # ---------------------
        nc.title = 'VLIDORT Simulation of GEOS-5 lidar background polarized reflectance'
        nc.institution = 'NASA/Goddard Space Flight Center'
        nc.source = 'Global Model and Assimilation Office'
        nc.history = 'VLIDORT simulation run on sampled GEOS-5'
        nc.satangles  = self.VZA
        nc.references = 'n/a'
        nc.contact = 'Patricia Castellanos <patricia.castellanos@nasa.gov>'
        nc.Conventions = 'CF'
        nc.inFile = self.inFile
        nc.comment = 'TOA reflectance simulations limited to SZA < 80'        
     
        # Create dimensions
        # -----------------
        nt = nc.createDimension('time',len(self.tyme))
        ls = nc.createDimension('ls',19)
        nz = nc.createDimension('lev',km)
        x  = nc.createDimension('x',1)
        y  = nc.createDimension('y',1)
        ch  = nc.createDimension('channel',1)

        # Coordinate variables
        # --------------------
        col = 'aer_Nv'
        if self.verbose: 
            print 'opening file',self.inFile.replace('%col',col)
        nctrj       = Dataset(self.inFile.replace('%col',col))        
        _copyVar(nctrj,nc,u'trjLon',dtype='f4',zlib=False,verbose=self.verbose)
        _copyVar(nctrj,nc,u'trjLat',dtype='f4',zlib=False,verbose=self.verbose)
        _copyVar(nctrj,nc,u'time', dtype='i4',zlib=False,verbose=self.verbose)
        _copyVar(nctrj,nc,u'lev', dtype='S1',zlib=False,verbose=self.verbose)
        _copyVar(nctrj,nc,u'isotime', dtype='S1',zlib=False,verbose=self.verbose)
        _copyVar(nctrj,nc,u'x',dtype='f4',zlib=False,verbose=self.verbose)
        _copyVar(nctrj,nc,u'y',dtype='f4',zlib=False,verbose=self.verbose)   
        nctrj.close()


        vza = nc.createVariable('sensor_zenith','f4',('time',),zlib=zlib)
        vza.long_name     = "sensor viewing zenith angle (VZA)"
        vza.missing_value = MISSING
        vza.units         = "degrees (positive forward view)"
        vza[:]            = self.VZA

        vaa = nc.createVariable('sensor_azimuth','f4',('time',),zlib=zlib)
        vaa.long_name     = "sensor viewing azimuth angle (VAA)"
        vaa.missing_value = MISSING
        vaa.units         = "degrees clockwise from North"
        
        vaa[:]            = self.VAA

        sza = nc.createVariable('solar_zenith','f4',('time',),zlib=zlib)
        sza.long_name     = "solar zenith angle (SZA)"
        sza.missing_value = MISSING
        sza.units         = "degrees"
        sza[:]            = self.SZA  

        saa = nc.createVariable('solar_azimuth','f4',('time',),zlib=zlib)
        saa.long_name     = "solar azimuth angle (SAA)"
        saa.missing_value = MISSING
        saa.units         = "degrees clockwise from North"
        saa[:]            = self.SAA


        # Write VLIDORT Outputs
        # ---------------------
        ref = nc.createVariable('toa_reflectance','f4',('time',),zlib=zlib,fill_value=MISSING)
        ref.standard_name = '%.2f nm TOA Reflectance' %self.channel
        ref.long_name     = '%.2f nm reflectance at the top of the atmosphere' %self.channel
        ref.missing_value = MISSING
        ref.units         = "None"
        ref[:]            = self.reflectance

        i = nc.createVariable('I','f4',('time',),zlib=zlib,fill_value=MISSING)
        i.standard_name = '%.2f nm TOA I' %self.channel
        i.long_name     = '%.2f nm sun normalized intensity at the top of the atmosphere' %self.channel
        i.missing_value = MISSING
        i.units         = "W m-2 sr-1 nm-1"
        i[:]            = self.I

        q = nc.createVariable('Q','f4',('time',),zlib=zlib,fill_value=MISSING)
        q.standard_name = '%.2f nm TOA Q' %self.channel
        q.long_name     = '%.2f nm Q-component of the stokes vector at the top of the atmopshere' %self.channel
        q.missing_value = MISSING
        q.units         = "W m-2 sr-1 nm-1"
        q[:]            = self.Q    

        u = nc.createVariable('U','f4',('time',),zlib=zlib,fill_value=MISSING)
        u.standard_name = '%.2f nm TOA U' %self.channel
        u.long_name     = '%.2f nm U-component of the stokes vector at the top of the atmopshere' %self.channel
        u.missing_value = MISSING
        u.units         = "W m-2 sr-1 nm-1"
        u[:]            = self.U

        sref = nc.createVariable('surf_reflectance','f4',('time',),zlib=zlib,fill_value=MISSING)
        sref.standard_name = '%.2f nm Surface Reflectance' %self.channel
        sref.long_name     = '%.2f nm Bi-Directional Surface Reflectance' %self.channel
        sref.missing_value = MISSING
        sref.units         = "None"
        sref[:]            = self.surf_reflectance

        sref = nc.createVariable('surf_reflectance_Q','f4',('time',),zlib=zlib,fill_value=MISSING)
        sref.standard_name = '%.2f nm Surface Reflectance Q' %self.channel
        sref.long_name     = '%.2f nm Bi-Directional Surface Reflectance Q' %self.channel
        sref.missing_value = MISSING
        sref.units         = "None"
        sref[:]            = self.BR_Q

        sref = nc.createVariable('surf_reflectance_U','f4',('time',),zlib=zlib,fill_value=MISSING)
        sref.standard_name = '%.2f nm Surface Reflectance U' %self.channel
        sref.long_name     = '%.2f nm Bi-Directional Surface Reflectance U' %self.channel
        sref.missing_value = MISSING
        sref.units         = "None"
        sref[:]            = self.BR_U

        rot = nc.createVariable('ROT','f4',('time','lev',),zlib=zlib,fill_value=MISSING)
        rot.long_name = '%.2f nm Rayleigh Optical Thickness' %self.channel
        rot.missing_value = MISSING
        rot.units         = "None"
        rot[:]            = self.ROT

        depol = nc.createVariable('rayleigh_depol_ratio','f4',('channel',),zlib=zlib,fill_value=MISSING)
        depol.long_name = '%.2f nm Rayleigh Depolrization Ratio' %self.channel
        depol.missing_value = MISSING
        depol.units         = "None"
        depol[:]            = self.depol_ratio

        mr = nc.createVariable('ocean_refractive_index','f4',('channel',),zlib=zlib,fill_value=MISSING)
        mr.long_name = '%.2f nm ocean refreactive index' %self.channel
        mr.missing_value = MISSING
        mr.units         = "None"
        mr[:]            = self.mr

        # Close the file
        # --------------
        nc.close()

        if self.verbose:
            print " <> wrote %s"%(self.outFile)

    
def get_chd(channel):
    chd = '%.2f'%channel
    chd = chd.replace('.','d')

    return chd
    
#------------------------------------ M A I N ------------------------------------

if __name__ == "__main__":
    date     = datetime(2006,01,01,00)
    nymd     = str(date.date()).replace('-','')
    hour     = str(date.hour).zfill(2)
    format   = 'NETCDF4_CLASSIC'

    rootDir  = '/discover/nobackup/projects/gmao/osse2/pub/c1440_NR/OBS/POLAR_LIDAR/CALIPSO/'

    inDir        = '{}/LevelB/Y{}/M{}'.format(rootDir,date.year,str(date.month).zfill(2))
    inFile       = '{}/calipso-g5nr.lb2.%col.{}_{}z.nc4'.format(inDir,nymd,hour)
    brdfDir      = '{}/LevelB/surface/BRDF/MCD43C1/006/Y{}/M{}'.format(rootDir,date.year,str(date.month).zfill(2))
    brdfFile     = '{}/calipso-g5nr.lb2.brdf.{}_{}z.nc4'.format(brdfDir,nymd,hour)
    ndviDir      = '{}/LevelB/surface/BPDF/NDVI/MYD13C2/006/Y{}/M{}'.format(rootDir,date.year,str(date.month).zfill(2))
    ndviFile     = '{}/calipso-g5nr.lb2.ndvi.{}_{}z.nc4'.format(ndviDir,nymd,hour)
    lcDir        = '{}/LevelB/surface/BPDF/LAND_COVER/MCD12C1/051/Y{}/M{}'.format(rootDir,date.year,str(date.month).zfill(2))
    lcFile       = '{}/calipso-g5nr.lb2.land_cover.{}_{}z.nc4'.format(lcDir,nymd,hour)    
    albedoType   = 'MODIS_BRDF_BPDF'
    VZAname      = 'POLDER'

    channel   = 550
    chd       = get_chd(channel)
    outDir    = '{}/LevelC/Y{}/M{}'.format(rootDir,date.year,str(date.month).zfill(2))
    outFile   = '{}/calipso-g5nr.vlidort.vector.{}.{}.{}_{}z_{}nm.nc4'.format(outDir,albedoType,VZAname,nymd,hour,chd)
    
    rcFile   = 'Aod_EOS.rc'
    orbit    = 'LEO'
    hgtss    = 705
    verbose  = True

    # Initialize VLIDORT class getting aerosol optical properties
    # -----------------------------------------------------------
    vlidort = LIDAR_VLIDORT(inFile,outFile,rcFile,
                            albedoType,
                            channel,                            
                            hgtss,
                            brdfFile=brdfFile,
                            ndviFile=ndviFile,
                            lcFile=lcFile,
                            verbose=verbose)

   
    # Run VLIDORT
    # if vlidort.nobs > 0:
    #     vlidort.runVLIDORT()
