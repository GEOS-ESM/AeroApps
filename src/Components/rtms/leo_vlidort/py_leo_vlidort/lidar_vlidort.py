#!/usr/bin/env python3

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

from datetime        import datetime, timedelta
from dateutil.parser import parse         as isoparser

from MAPL.constants import *
import LidarAngles_    
import VLIDORT_POLAR_ 
from .copyvar  import _copyVar
from scipy import interpolate
from   MAPL  import config

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

WrapperFuncs = {'MODIS_BRDF'     : VLIDORT_POLAR_.vector_brdf_modis,
                'MODIS_BRDF_BPDF': VLIDORT_POLAR_.vector_brdf_modis_bpdf,
                'LAMBERTIAN'     : VLIDORT_POLAR_.vector_lambert,
                'GissCX'         : VLIDORT_POLAR_.vector_gisscx,
                'CX'             : VLIDORT_POLAR_.vector_cx,
                'ROT_CALC'       : VLIDORT_POLAR_.rot_calc}   


MISSING = -1.e+20

class LIDAR_VLIDORT(object):
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
        albedoReader = getattr(self,SurfaceFuncs['CX'])
        albedoReader()

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
    def BPDFinputs(self):
        """
        Read in NDVI and Landuse Coefficient- should have been sampled already
        """
        if self.verbose:
            print('opening file',self.ndviFile)

        nc = Dataset(self.ndviFile)
        NDVI = nc.variables['NDVI'][:]
        I = NDVI < -900
        NDVI[I] = MISSING
        nc.close()

        if self.verbose:
            print('opening file',self.lcFile)
        nc = Dataset(self.lcFile)
        BPDFcoef = nc.variables['BPDFcoef'][:]
        I = BPDFcoef < -900
        BPDFcoef[I] = MISSING
        nc.close()

        self.iGoodLand = self.iGoodLand & (NDVI != MISSING)
        self.nobsLand  = np.sum(self.iGoodLand)

        #BPDFparam(nparam,nch,nobs)
        self.BPDFparam = np.zeros([3,1,len(self.tyme)])
        self.BPDFparam[0,0,:] = 1.5
        self.BPDFparam[1,0,:] = NDVI
        self.BPDFparam[2,0,:] = BPDFcoef        

    # ---
    def LandSeaMask(self):
        """
        Read in invariant dataset
        """
        col = 'asm_Nx'
        if self.verbose: 
            print('opening file',self.inFile.replace('%col',col))
        nc       = Dataset(self.inFile.replace('%col',col))

        for sds in self.SDS_INV:
            sds_ = sds
            if sds in ncALIAS:
                sds_ = ncALIAS[sds]
            var = nc.variables[sds_][:]
            self.__dict__[sds].append(var)   

        # Make lists into arrays
        for sds in self.SDS_INV:
            self.__dict__[sds] = np.concatenate(self.__dict__[sds])            

        self.iLand = (self.FRLAND >= 0.99) & self.iGoodLand
        self.iSea  = self.FRLAND < 0.99

        # self.iGood = self.iGood & iGood
        self.nobsLand  = np.sum(self.iGood & self.iLand)
        self.nobsSea   = np.sum(self.iGood & self.iSea)

    # --
    def computeAtmos(self):

        pe, ze, te = getEdgeVars(self)

        self.pe = pe # (km,nobs)
        self.ze = ze 
        self.te = te

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
            var = nc.variables[sds_][:]
            self.__dict__[sds].append(var)

        if len(self.SDS_MET) > 0:
            col = 'met_Nv'
            if self.verbose: 
                print('opening file',self.inFile.replace('%col',col))        
            nc       = Dataset(self.inFile.replace('%col',col))

            for sds in self.SDS_MET:
                sds_ = sds
                if sds in ncALIAS:
                    sds_ = ncALIAS[sds]
                var = nc.variables[sds_][:]
                self.__dict__[sds].append(var)

    def readSampledWindCX(self):
        col = 'met_Nv'
        if self.verbose: 
            print('opening file',self.inFile.replace('%col',col))        
        nc       = Dataset(self.inFile.replace('%col',col))

        for sds in self.SDS_CX:
            sds_ = sds
            if sds in ncALIAS:
                sds_ = ncALIAS[sds]
            var = nc.variables[sds_][:]
            self.__dict__[sds].append(var)

        for sds in self.SDS_CX:
            self.__dict__[sds] = np.concatenate(self.__dict__[sds])

        self.mr = np.array([1.333])


    # --- 
    def readSampledMODISBRDF(self):
        """
        Read in MODIS BRDF kernel weights
        that have already been sampled on lidar track
        """
        MODIS_channels = np.array(['470','550','650','850','1200','1600','2100'])
        chs = str(int(self.channel))
        if chs in MODIS_channels:
            SDS = 'Riso'+chs,'Rgeo'+chs,'Rvol'+chs
        else:
            dch = self.channel - np.array(MODIS_channels).astype('int')
            chmin = np.argmin(dch[dch>0])
            chmin = MODIS_channels[dch>0][chmin]
            chmax = np.argmax(dch[dch<0])
            chmax = MODIS_channels[dch<0][chmax]

            SDS = 'Riso'+chmin,'Rgeo'+chmin,'Rvol'+chmin,'Riso'+chmax,'Rgeo'+chmax,'Rvol'+chmax

        if self.verbose:
            print('opening BRDF file ',self.brdfFile)
        nc = Dataset(self.brdfFile)

        for sds in SDS:
            self.__dict__[sds] = np.array(nc.variables[sds][:])

        missing_value = nc.variables[sds].missing_value
        nc.close()

        nobs = len(self.__dict__[sds])
        # Interpolate if necessary
        if chs not in MODIS_channels:
            X = np.array([chmin,chmax]).astype('int')
            for R in ['Riso','Rgeo','Rvol']:
                sds = R + chs
                self.__dict__[sds] = np.empty([nobs])
                for i in range(nobs):
                    Y = np.array([self.__dict__[R+chmin][i],self.__dict__[R+chmax][i]])
                    if missing_value in Y:
                        self.__dict__[sds][i] = missing_value
                    else:
                        f = interpolate.interp1d(X, Y)
                        self.__dict__[sds][i] = f([int(chs)])

            SDS = 'Riso'+chs,'Rgeo'+chs,'Rvol'+chs
        
        # Check for missing kernel weights
        Riso = self.__dict__['Riso'+chs]
        Rgeo = self.__dict__['Rgeo'+chs]
        Rvol = self.__dict__['Rvol'+chs]
        iGood = (Riso != missing_value) & (Rgeo != missing_value) & (Rvol != missing_value)
        self.iGoodLand = self.iGoodLand & iGood
        self.nobsLand = np.sum(self.iGoodLand)

        for sds in SDS:
            self.__dict__[sds].shape = (1,1,nobs)

        # [nkernel,nch,nobs]
        self.kernel_wt = np.append(self.__dict__['Riso'+chs],self.__dict__['Rgeo'+chs],axis=0)
        self.kernel_wt = np.append(self.kernel_wt,self.__dict__['Rvol'+chs],axis=0)

        param1 = np.array([2]*nobs)
        param2 = np.array([1]*nobs)

        param1.shape = (1,1,nobs)
        param2.shape = (1,1,nobs)

        # [nparam,nch,nobs]
        self.RTLSparam = np.append(param1,param2,axis=0)

    # --- 
    def readHybridMODISBRDF(self):
        """
        Read in MODIS BRDF kernel weights
        that have already been sampled on lidar track
        for wavelngths between 388 and 470
        """
        MODIS_channels = np.array(['470','550','650','850','1200','1600','2100'])
        LER_channels   = np.array(['354','388'])
        MODISchs = '470' 
        mSDS = 'Riso'+MODISchs,'Rgeo'+MODISchs,'Rvol'+MODISchs
        LERchs = '388'
        lSDS = 'SRFLER' + LERchs
        chs = str(int(self.channel))

        if self.verbose:
            print('opening BRDF abledo file ',self.brdfFile)
        nc = Dataset(self.brdfFile)

        for sds in mSDS:
            self.__dict__[sds] = nc.variables[sds][:]

        missing_value = nc.variables[sds].missing_value
        nc.close()

        if self.verbose:
            print('opening LER albedo file ',self.lerFile)
        nc = Dataset(self.lerFile)

        self.__dict__[lSDS] = nc.variables[lSDS][:]
        nc.close()

        nobs = len(self.__dict__[sds])
        # Interpolate 
        X = np.array([388,470])
        #Riso
        sds = 'Riso' + chs
        R   = 'Riso' + MODISchs
        self.__dict__[sds] = np.empty([nobs])
        for i in range(nobs):
            Y = np.array([self.__dict__[lSDS][i],self.__dict__[R][i]])
            if missing_value in Y:
                self.__dict__[sds][i] = missing_value
            else:
                f = interpolate.interp1d(X, Y)
                self.__dict__[sds][i] = f([int(chs)])

        #Rgeo and Rvol
        for R in 'Rgeo','Rvol':
            sds = R + chs
            self.__dict__[sds] = np.empty([nobs])
            for i in range(nobs):
                Y = np.array([0.0,self.__dict__[R+MODISchs][i]])
                if missing_value in Y:
                    self.__dict__[sds][i] = missing_value
                else:
                    f = interpolate.interp1d(X, Y)
                    self.__dict__[sds][i] = f([int(chs)])

        SDS = 'Riso'+chs,'Rgeo'+chs,'Rvol'+chs
        
        # Check for missing kernel weights
        Riso = self.__dict__['Riso'+chs]
        Rgeo = self.__dict__['Rgeo'+chs]
        Rvol = self.__dict__['Rvol'+chs]
        iGood = (Riso != missing_value) & (Rgeo != missing_value) & (Rvol != missing_value)
        self.iGoodLand = self.iGoodLand & iGood
        self.nobsLand = np.sum(self.iGoodLand)

        for sds in SDS:
            self.__dict__[sds].shape = (1,1,nobs)

        # [nkernel,nch,nobs]
        self.kernel_wt = np.append(self.__dict__['Riso'+chs],self.__dict__['Rgeo'+chs],axis=0)
        self.kernel_wt = np.append(self.kernel_wt,self.__dict__['Rvol'+chs],axis=0)

        param1 = np.array([2]*nobs)
        param2 = np.array([1]*nobs)

        param1.shape = (1,1,nobs)
        param2.shape = (1,1,nobs)

        # [nparam,nch,nobs]
        self.RTLSparam = np.append(param1,param2,axis=0)

    def readSampledLER(self):
        """
        Read in sampler LER files.
        Fill albedo attribute
        """
        LER_channels   = np.array(['354','388'])   
        chs = str(int(self.channel))

        SDS = []
        for lchs in LER_channels:
            SDS.append('SRFLER'+lchs)

        if self.verbose:
            print('opening LER albedo file ',self.lerFile)
        nc = Dataset(self.lerFile)

        for sds in SDS:
            self.__dict__[sds] = np.squeeze(nc.variables[sds][:])

        missing_value = nc.variables[sds].missing_value
        nobs = len(self.__dict__[sds])
        nc.close()

        sds = 'SRFLER'+chs
        # interpolate if needed
        if chs not in LER_channels:
            dch   = self.channel - np.array(LER_channels).astype('int')
            chmin = np.argmin(dch[dch>0])
            chmin = LER_channels[dch>0][chmin]
            chmax = np.argmax(dch[dch<0])
            chmax = LER_channels[dch<0][chmax]            
            X = np.array([chmin,chmax]).astype('int')

            self.__dict__[sds] = np.empty([nobs])
            for i in range(nobs):
                Y = np.array([self.__dict__['SRFLER'+chmin][i],self.__dict__['SRFLER'+chmax][i]])
                if missing_value in Y:
                    self.__dict__[sds][i] = missing_value
                else:
                    f = interpolate.interp1d(X, Y)
                    self.__dict__[sds][i] = f([int(chs)]) 


        # Check for missing values
        iGood = (self.__dict__[sds] != missing_value) 
        self.iGoodLand= self.iGoodLand & iGood
        self.nobsLand = np.sum(self.iGoodLand)

        # (nobs,nch)
        self.__dict__[sds].shape = (nobs,1) 


        self.albedo = self.__dict__[sds]     

    #---
    def computeMie(self):
        """
        Computes aerosol optical quantities 
        """
        tau,ssa,g,pmom = getAOPvector(self,self.channel,
                                 vnames=self.AERNAMES,
                                 Verbose=True,
                                 rcfile=self.rcFile,
                                 nMom=self.nMom)
        self.tau = tau  #(km,nch,nobs)
        self.ssa = ssa  #(km,nch,nobs)
        self.g   = g    #(km,nch,nobs)
        self.pmom = pmom  #(km,nch,nobs,nMom,nPol)
        # Multiply by -1 to go from Mischenko convention to VLIDORT
        pmom[:,:,:,:,1] = -1.*pmom[:,:,:,:,1]
        pmom[:,:,:,:,3] = -1.*pmom[:,:,:,:,3]


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

        
        # loop though LAND and SEA
        for surface in surfList:

            iGood = self.__dict__['i'+surface]
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

            # Get VLIDORT wrapper name from dictionary
            vlidortWrapper = WrapperFuncs[albedoType]

            # run VLIDORT
            runFunction = self.__getattribute__(albedoType+'_run')
            I,Q,U,reflectance,surf_reflectance,BR_Q,BR_U = runFunction(vlidortWrapper,rot,depol_ratio,tau,ssa,pmom,pe,ze,te,sza,raa,vza,iGood)
                                                
            self.I[iGood] = np.squeeze(I)
            self.reflectance[iGood] = np.squeeze(reflectance)
            self.surf_reflectance[iGood] = np.squeeze(surf_reflectance)
            self.Q[iGood] = np.squeeze(Q)
            self.U[iGood] = np.squeeze(U) 
            self.BR_Q[iGood] = np.squeeze(BR_Q)
            self.BR_U[iGood] = np.squeeze(BR_U) 


        self.writeNC()
    #---
    def CX_run(self,vlidortWrapper,ROT,depol_ratio,tau,ssa,pmom,pe,ze,te,sza,raa,vza,iGood):
        u10m = self.U10M[iGood]
        v10m = self.V10M[iGood]

        args = [self.channel, self.nstreams, self.plane_parallel, ROT, depol_ratio, tau, ssa, pmom, 
                pe, ze, te, 
                u10m, v10m, self.mr, 
                sza, raa, vza, 
                MISSING,
                self.verbose]

        # Call VLIDORT wrapper function
        I, reflectance, surf_reflectance, Q, U, BR_Q, BR_U, rc = vlidortWrapper(*args)                        

        return I,Q,U,reflectance,surf_reflectance,BR_Q,BR_U
    #---
    def MODIS_BRDF_run(self,vlidortWrapper,ROT,depol_ratio,tau,ssa,pmom,pe,ze,te,sza,raa,vza,iGood):
        kernel_wt = self.kernel_wt[:,:,iGood]
        param     = self.RTLSparam[:,:,iGood]                
        

        args = [self.channel, self.nstreams, self.plane_parallel, ROT, depol_ratio, tau, ssa, pmom, 
                pe, ze, te, 
                kernel_wt, param, 
                sza, raa, vza, 
                MISSING,
                self.verbose]

        # Call VLIDORT wrapper function
        I, reflectance, surf_reflectance, Q, U, BR_Q, BR_U, rc = vlidortWrapper(*args)                        

        return I,Q,U,reflectance,surf_reflectance,BR_Q,BR_U       
    #---    
    def MODIS_BRDF_BPDF_run(self,vlidortWrapper,ROT,depol_ratio,tau,ssa,pmom,pe,ze,te,sza,raa,vza,iGood):
        # For albedo
        kernel_wt = self.kernel_wt[:,:,iGood]
        RTLSparam = self.RTLSparam[:,:,iGood] 
        RTLSparam = np.append(RTLSparam,np.zeros([1,1,self.nobsLand]),axis=0) 

        # For BPDF
        BPDFparam = self.BPDFparam[:,:,iGood]

        # Loop through one by one
        # Some land covers do not have polarization (i.e. urban)
        I                = np.zeros([self.nobsLand,1])
        Q                = np.zeros([self.nobsLand,1])
        U                = np.zeros([self.nobsLand,1])
        reflectance      = np.zeros([self.nobsLand,1])
        surf_reflectance = np.zeros([self.nobsLand,1])
        BR_Q             = np.zeros([self.nobsLand,1])
        BR_U             = np.zeros([self.nobsLand,1])
        nlev             = tau.shape[0]

        for p in range(self.nobsLand):
            if BPDFparam[2,0,p] == MISSING: 
                args = [self.channel,
                        self.nstreams,
                        self.plane_parallel,
                        ROT[:,p:p+1,:],
                        depol_ratio,
                        tau[:,:,p:p+1], 
                        ssa[:,:,p:p+1], 
                        pmom[:,:,p:p+1,:,:], 
                        pe[:,p:p+1], 
                        ze[:,p:p+1], 
                        te[:,p:p+1], 
                        kernel_wt[:,:,p:p+1], 
                        RTLSparam[:,:,p:p+1], 
                        sza[p:p+1], 
                        raa[p:p+1], 
                        vza[p:p+1], 
                        MISSING,
                        self.verbose]

                BRDFvlidortWrapper = WrapperFuncs['MODIS_BRDF']
                # Call VLIDORT wrapper function
                I_, reflectance_, surf_reflectance_, Q_, U_, BR_Q_, BR_U_, rc = BRDFvlidortWrapper(*args)                         

            else:
                args = [self.channel,
                        self.nstreams,
                        self.plane_parallel,
                        ROT[:,p:p+1,:],
                        depol_ratio,
                        tau[:,:,p:p+1], 
                        ssa[:,:,p:p+1], 
                        pmom[:,:,p:p+1,:,:], 
                        pe[:,p:p+1], 
                        ze[:,p:p+1], 
                        te[:,p:p+1], 
                        kernel_wt[:,:,p:p+1], 
                        RTLSparam[:,:,p:p+1], 
                        BPDFparam[:,:,p:p+1],
                        sza[p:p+1], 
                        raa[p:p+1], 
                        vza[p:p+1], 
                        MISSING,
                        self.verbose]

                # Call VLIDORT wrapper function
                I_, reflectance_, surf_reflectance_, Q_, U_, BR_Q_, BR_U_, rc = vlidortWrapper(*args)
    
            I[p:p+1,:] = I_
            Q[p:p+1,:] = Q_
            U[p:p+1,:] = U_
            reflectance[p:p+1,:] = reflectance_
            surf_reflectance[p:p+1,:] = surf_reflectance_
            BR_Q[p:p+1,:] = BR_Q_
            BR_U[p:p+1,:] = BR_U_

        return I,Q,U,reflectance,surf_reflectance,BR_Q,BR_U
    #---
    def LAMBERTIAN_run(self,vlidortWrapper,ROT,depol_ratio,tau,ssa,pmom,pe,ze,te,sza,raa,vza,iGood):
        albedo = self.albedo[iGood,:]
        
        args = [self.channel,self.nstreams,self.plane_parallel,ROT,depol_ratio,tau, ssa, pmom, 
                pe, ze, te, 
                albedo, 
                sza, raa, vza, 
                MISSING,
                self.verbose]

        # Call VLIDORT wrapper function
        I, reflectance, Q, U, rc = vlidortWrapper(*args)  
        surf_reflectance = albedo            

        BR_Q = None
        BR_U = None
        return I,Q,U,reflectance,surf_reflectance,BR_Q,BR_U
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
            print('opening file',self.inFile.replace('%col',col))
        nctrj       = Dataset(self.inFile.replace('%col',col))        
        _copyVar(nctrj,nc,'trjLon',dtype='f4',zlib=False,verbose=self.verbose)
        _copyVar(nctrj,nc,'trjLat',dtype='f4',zlib=False,verbose=self.verbose)
        _copyVar(nctrj,nc,'time', dtype='i4',zlib=False,verbose=self.verbose)
        _copyVar(nctrj,nc,'lev', dtype='S1',zlib=False,verbose=self.verbose)
        _copyVar(nctrj,nc,'isotime', dtype='S1',zlib=False,verbose=self.verbose)
        _copyVar(nctrj,nc,'x',dtype='f4',zlib=False,verbose=self.verbose)
        _copyVar(nctrj,nc,'y',dtype='f4',zlib=False,verbose=self.verbose)   
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
            print(" <> wrote %s"%(self.outFile))

    
def get_chd(channel):
    chd = '%.2f'%channel
    chd = chd.replace('.','d')

    return chd
    
#------------------------------------ M A I N ------------------------------------

if __name__ == "__main__":
    date     = datetime(2006,0o1,0o1,00)
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

    channel   = 532
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
