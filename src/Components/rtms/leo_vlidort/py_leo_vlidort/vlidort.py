#!/usr/bin/env python3

"""
    Parent class with utilities needed to 
    setup input data
    and run python wrappers for vlidort
"""

import os
from   netCDF4 import Dataset
from   mieobs  import  getAOPvector, getEdgeVars
from   leo_vlidort import VLIDORT_POLAR_ 
import numpy   as np
from MAPL.constants import *
from scipy import interpolate

MISSING = -1.e+20

ncALIAS = {'LONGITUDE': 'trjLon',
           'LATITUDE': 'trjLat'}

WrapperFuncs = {'MODIS_BRDF'     : VLIDORT_POLAR_.vector_brdf_modis,
                'MODIS_BRDF_CLOUD'     : VLIDORT_POLAR_.vector_brdf_modis_cloud,
                'MODIS_BRDF_BPDF': VLIDORT_POLAR_.vector_brdf_modis_bpdf,
                'LAMBERTIAN'     : VLIDORT_POLAR_.vector_lambert,
                'LAMBERTIAN_CLOUD'     : VLIDORT_POLAR_.vector_lambert_cloud,
                'GissCX'         : VLIDORT_POLAR_.vector_gisscx,
                'CX'             : VLIDORT_POLAR_.vector_cx,
                'CX_CLOUD'             : VLIDORT_POLAR_.vector_cx_cloud,
                'ROT_CALC'       : VLIDORT_POLAR_.rot_calc}              

#---
def CX_run(args):
       
    # Call VLIDORT wrapper function
    I, reflectance, surf_reflectance, Q, U, BR_Q, BR_U, rc = VLIDORT_POLAR_.vector_cx(*args)

    return I,Q,U,reflectance,surf_reflectance,BR_Q,BR_U
#---
def CX_CLOUD_run(args):

    # Call VLIDORT wrapper function
    I, reflectance, surf_reflectance, Q, U, BR_Q, BR_U, rc = VLIDORT_POLAR_.vector_cx_cloud(*args)

    return I,Q,U,reflectance,surf_reflectance,BR_Q,BR_U

#---
def LAMBERTIAN_run(args):

    # Call VLIDORT wrapper function
    I, reflectance, Q, U, rc = VLIDORT_POLAR_.vector_lambert(*args)

    surf_reflectance = None
    BR_Q = None
    BR_U = None
    return I,Q,U,reflectance,surf_reflectance,BR_Q,BR_U
#---
def LAMBERTIAN_CLOUD_run(args):

    # Call VLIDORT wrapper function
    I, reflectance, Q, U, rc = VLIDORT_POLAR_.vector_lambert_cloud(*args)

    surf_reflectance = None
    BR_Q = None
    BR_U = None
    return I,Q,U,reflectance,surf_reflectance,BR_Q,BR_U

#---
def MODIS_BRDF_run(args):

    # Call VLIDORT wrapper function
    I, reflectance, surf_reflectance, Q, U, BR_Q, BR_U, rc = VLIDORT_POLAR_.vector_brdf_modis(*args)

    return I,Q,U,reflectance,surf_reflectance,BR_Q,BR_U
#---
def MODIS_BRDF_CLOUD_run(args):

    # Call VLIDORT wrapper function
    I, reflectance, surf_reflectance, Q, U, BR_Q, BR_U, rc = VLIDORT_POLAR_.vector_brdf_modis_cloud(*args)

    return I,Q,U,reflectance,surf_reflectance,BR_Q,BR_U

class VLIDORT(object):
    """ 
    Utilities for setting up input data 
    and running python wrapper for vlidort
    """
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
            # cast from masked array
            self.__dict__[sds] = np.array(nc.variables[sds][:])

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

    # ---
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
    def readIRR(self):
        """ 
        Read in solar irradiance
        """
        nch = 1
        self.flux_factor = np.ones([nch,self.nobs])
     
    #---
    def computeAlpha(self):
        """
        Computes trace gas absorption optical depth
        """
        km = 72
        nch = 1
        self.alpha = np.zeros([km,self.nobs,nch])

    #---
    def computeMie(self):
        """
        Computes aerosol optical quantities 
        """
        tau,ssa,g,pmom = getAOPvector(self,self.channel,
                                 vnames=self.AERNAMES,
                                 Verbose=self.verbose,
                                 rcfile=self.rcFile,
                                 nMom=self.nMom)
        self.tau = tau  #(km,nch,nobs)
        self.ssa = ssa  #(km,nch,nobs)
        self.g   = g    #(km,nch,nobs)

        # Multiply by -1 to go from Mischenko convention to VLIDORT
        pmom[:,:,:,:,1] = -1.*pmom[:,:,:,:,1]
        pmom[:,:,:,:,3] = -1.*pmom[:,:,:,:,3]
        self.pmom = pmom #(km,nch,nobs,nMom,nPol)
    # --
    def computeAtmos(self):

        pe, ze, te = getEdgeVars(self)

        self.pe = pe # (km,nobs)
        self.ze = ze 
        self.te = te

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

        BR_Q = 0.
        BR_U = 0.
        return I,Q,U,reflectance,surf_reflectance,BR_Q,BR_U
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
        

def get_chd(channel):
    chd = '%.2f'%channel
    chd = chd.replace('.','d')

    return chd
