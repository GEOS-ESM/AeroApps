#!/usr/bin/env python

"""
    Calculates polarized TOA radiance for a multiangle polarimeter on a lidar track.
    Model fields have already been sampled using trj_sampler

    Adapted from ext_sampler.py
    Patricia Castellanos, May, 2017

"""

import os
import MieObs_
from   netCDF4 import Dataset
from   mieobs  import  getAOPvector, getEdgeVars
import numpy   as np

from datetime        import datetime, timedelta
from dateutil.parser import parse         as isoparser

from MAPL.constants import *
import LidarAngles_    
import VLIDORT_ 
from copyvar  import _copyVar
from scipy import interpolate

# Generic Lists of Varnames and Units
VNAMES_DU = ['DU001','DU002','DU003','DU004','DU005']
VNAMES_SS = ['SS001','SS002','SS003','SS004','SS005']
VNAMES_BC = ['BCPHOBIC','BCPHILIC']
VNAMES_OC = ['OCPHOBIC','OCPHILIC']
VNAMES_SU = ['SO4']
MieVarsNames = ['ext','scatext','backscat','aback_sfc','aback_toa','depol','ext2back','tau','ssa','g']
MieVarsUnits = ['km-1','km-1','km-1 sr-1','sr-1','sr-1','unitless','sr','unitless','unitless','unitless']

META    = ['DELP','PS','RH','AIRDENS','LONGITUDE','LATITUDE','isotime']
AERNAMES = VNAMES_SU + VNAMES_SS + VNAMES_OC + VNAMES_BC + VNAMES_DU
SDS_AER = META + AERNAMES
SDS_MET = []#['CLDTOT']
SDS_INV = ['FRLAND']

ncALIAS = {'LONGITUDE': 'trjLon',
           'LATITUDE': 'trjLat'}


nMom     = 300

VZAdic = {'POLDER': np.array([3.66, 11., 18.33, 25.66, 33, 40.33, 47.66, 55.0])}

HGTdic = {'LEO': 705,
          'ISS': 400}


SurfaceFuncs = {'MODIS_BRDF'     : 'readSampledMODISBRDF',
                'MODIS_BRDF_BPDF': 'readSampledMODISBRDF',
                'LAMBERTIAN'     : 'readSampledLER'}

WrapperFuncs = {'MODIS_BRDF'     : VLIDORT_.vector_brdf_modis,
                'MODIS_BRDF_BPDF': VLIDORT_.vector_brdf_modis_bpdf,
                'LAMBERTIAN'     : VLIDORT_.vector_lambert}   

LandAlbedos  = 'MODIS_BRDF','MODIS_BRDF_BPDF','LAMBERTIAN'


MISSING = -1.e+20

class POLAR_VLIDORT(object):
    """
    Everything needed for calling VLIDORT
    GEOS-5 has already been sampled on lidar track
    """
    def __init__(self,inFile,outFile,rcFile,albedoFile,albedoType,
                channel,VZA,hgtss,
                ndviFile=None,
                lcFile=None,
                lerFile=None,
                verbose=False,
                extOnly=False):
        self.SDS_AER = SDS_AER
        self.SDS_MET = SDS_MET
        self.SDS_INV = SDS_INV
        self.AERNAMES = AERNAMES
        self.inFile  = inFile
        self.outFile = outFile
        self.albedoFile = albedoFile
        self.albedoType = albedoType
        self.rcFile  = rcFile
        self.channel = channel
        self.verbose = verbose
        self.nMom    = nMom
        self.lcFile = lcFile
        self.ndviFile = ndviFile
        self.lerFile  = lerFile

        # initialize empty lists
        for sds in self.SDS_AER+self.SDS_MET+self.SDS_INV:
            self.__dict__[sds] = []

        # Read in data model data
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

        if not extOnly:
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


            # Calculate aerosol optical properties
            self.computeMie()

            # Calculate atmospheric profile properties needed for Rayleigh calc
            self.computeAtmos()

        # Calculate Scene Geometry
        self.VZA = VZA
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
            print 'opening file',self.ndviFile

        nc = Dataset(self.ndviFile)
        NDVI = nc.variables['NDVI'][:]
        I = NDVI < -900
        NDVI[I] = MISSING
        nc.close()

        if self.verbose:
            print 'opening file',self.lcFile
        nc = Dataset(self.lcFile)
        BPDFcoef = nc.variables['BPDFcoef'][:]
        I = BPDFcoef < -900
        BPDFcoef[I] = MISSING
        nc.close()

        self.iGood = self.iGood & (NDVI != MISSING)
        self.nobs  = np.sum(self.iGood)

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
            print 'opening file',self.inFile.replace('%col',col)
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

        if self.albedoType in LandAlbedos:
            iGood = self.FRLAND >= 0.99
        else:
            iGood = self.FRLAND < 0.99

        self.iGood = self.iGood & iGood
        self.nobs  = np.sum(self.iGood)

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
        VAAf  = []
        VAAb  = []
        for i,tyme in enumerate(self.tyme):
            if i == len(self.tyme)-1:
                CLAT = self.LATITUDE[i-1]
                CLON = self.LONGITUDE[i-1]
            else:
                CLAT = self.LATITUDE[i+1]
                CLON = self.LONGITUDE[i+1]
            
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

            if i == len(self.tyme)-1:
                VAAb.append(sat_angles[0][0])
                vaaf = sat_angles[0][0] - 180.0
                if vaaf < 0:
                    vaaf = vaaf + 360.
                VAAf.append(vaaf)

            else:
                VAAf.append(sat_angles[0][0])
                vaab = sat_angles[0][0]+180
                if vaab > 360.:
                    vaab = vaab - 360.
                VAAb.append(vaab)


        self.SZA = np.array(SZA)
        self.SAA = np.array(SAA)
        self.VAAf = np.array(VAAf)
        self.VAAb = np.array(VAAb)

        # define RAA according to photon travel direction
        saa = self.SAA + 180.0
        I = saa >= 360.
        saa[I] = saa[I] - 360.

        self.RAAb = self.VAAb - saa
        self.RAAf = self.VAAf - saa

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
            print 'opening file',self.inFile.replace('%col',col)
        nc       = Dataset(self.inFile.replace('%col',col))

        for sds in self.SDS_AER:
            sds_ = sds
            if sds in ncALIAS:
                sds_ = ncALIAS[sds]
            var = nc.variables[sds_][:]
            self.__dict__[sds].append(var)

        # col = 'met_Nv'
        # if self.verbose: 
        #     print 'opening file',self.inFile.replace('%col',col)        
        # nc       = Dataset(self.inFile.replace('%col',col))

        # for sds in self.SDS_MET:
        #     sds_ = sds
        #     if sds in ncALIAS:
        #         sds_ = ncALIAS[sds]
        #     var = nc.variables[sds_][:]
        #     self.__dict__[sds].append(var)

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
            chmin = np.argmax(dch[dch<0])
            chmin = MODIS_channels[dch<0][chmin]
            chmax = np.argmin(dch[dch>0])
            chmax = MODIS_channels[dch>0][chmax]

            SDS = 'Riso'+chmin,'Rgeo'+chmin,'Rvol'+chmin,'Riso'+chmax,'Rgeo'+chmax,'Rvol'+chmax

        if self.verbose:
            print 'opening abledo file ',self.albedoFile
        nc = Dataset(self.albedoFile)

        for sds in SDS:
            self.__dict__[sds] = nc.variables[sds][:]

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
        self.iGood = self.iGood & iGood
        self.nobs = np.sum(self.iGood)

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
            print 'opening BRDF abledo file ',self.albedoFile
        nc = Dataset(self.albedoFile)

        for sds in mSDS:
            self.__dict__[sds] = nc.variables[sds][:]

        missing_value = nc.variables[sds].missing_value
        nc.close()

        if self.verbose:
            print 'opening LER albedo file ',self.lerFile
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
        self.iGood = self.iGood & iGood
        self.nobs = np.sum(self.iGood)

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
            print 'opening LER albedo file ',self.lerFile
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
            chmin = np.argmax(dch[dch<0])
            chmin = LER_channels[dch<0][chmin]
            chmax = np.argmin(dch[dch>0])
            chmax = LER_channels[dch>0][chmax]            
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
        self.iGood = self.iGood & iGood
        self.nobs = np.sum(self.iGood)

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

    # --
    def runExt(self):
        """
        run ext_sampler.py 
        """
        col = 'aer_Nv'
        if self.verbose: 
            print 'running ext_sampler on file',self.inFile.replace('%col',col)

        outDir = os.path.dirname(self.outFile)
        instname = os.path.basename(self.inFile).split('.')[0]
        date_ch   = os.path.basename(self.inFile).split('.')[-2]
        outFile = '{}/{}.lc2.ext.{}_{}nm.nc'.format(outDir,instname,date_ch,get_chd(self.channel))
        Options =     " --input=" + self.inFile.replace('%col',col)      + \
                      " --output=" + outFile       + \
                      " --rc=" + self.rcFile      + \
                      " --format=NETCDF4_CLASSIC"      + \
                      " --channel=%d" %self.channel + \
                      " --intensive"     
                      

        if not os.path.exists(os.path.dirname(outFile)):
            os.makedirs(os.path.dirname(outFile))

        cmd = 'ext_sampler.py {} '.format(Options)  

        if os.system(cmd):
            raise ValueError, "ext_sampler.py failed for %s "%(self.inFile.replace('%col',col))       


    def runVLIDORT(self):
        """
        Calls VLIDORT 
        """

        # Only do good obs
        # self.iGood = np.arange(len(self.iGood))[self.iGood][0:1]
        # self.nobs  = 1
        sza  = self.SZA[self.iGood]
        raaf = self.RAAf[self.iGood]
        raab = self.RAAb[self.iGood]

        tau = self.tau[:,:,self.iGood]
        ssa = self.ssa[:,:,self.iGood]
        pmom = self.pmom[:,:,self.iGood,:,:]
        pe   = self.pe[:,self.iGood]
        ze   = self.ze[:,self.iGood]
        te   = self.te[:,self.iGood]

        # Initiate output arrays
        nangles = len(self.VZA)
        ntime   = len(self.tyme)
        nlev    = tau.shape[0]
        self.I = np.ones([ntime,2*nangles])*MISSING
        self.Q = np.ones([ntime,2*nangles])*MISSING
        self.U = np.ones([ntime,2*nangles])*MISSING
        self.reflectance = np.ones([ntime,2*nangles])*MISSING
        self.surf_reflectance = np.ones([ntime,2*nangles])*MISSING
        self.ROT = np.ones([ntime,nlev])*MISSING

        #backward directions, forward directions
        VZA = np.append(self.VZA[::-1],self.VZA)

        raaf.shape = raaf.shape + (1,)
        raaf = np.repeat(raaf,len(self.VZA),axis=1)
        
        raab.shape = raab.shape + (1,)        
        raab = np.repeat(raab,len(self.VZA),axis=1)

        RAA = np.append(raab,raaf,axis=1)
        
        # Get VLIDORT wrapper name from dictionary
        vlidortWrapper = WrapperFuncs[self.albedoType]

        # loop through viewing angles and run VLIDORT
        for i in range(2*nangles):
            vza  = np.array([VZA[i]]*self.nobs)
            raa  = RAA[:,i]
            
            # get args list for each surface model
            if self.albedoType == 'MODIS_BRDF':
                kernel_wt = self.kernel_wt[:,:,self.iGood]
                param     = self.RTLSparam[:,:,self.iGood]                
                

                args = [self.channel,tau, ssa, pmom, 
                        pe, ze, te, 
                        kernel_wt, param, 
                        sza, raa, vza, 
                        MISSING,
                        self.verbose]

                # Call VLIDORT wrapper function
                I, reflectance, ROT, surf_reflectance, Q, U, rc = vlidortWrapper(*args)                        
                
            elif self.albedoType == 'MODIS_BRDF_BPDF':
                # For albedo
                kernel_wt = self.kernel_wt[:,:,self.iGood]
                RTLSparam = self.RTLSparam[:,:,self.iGood] 
                RTLSparam = np.append(RTLSparam,np.zeros([1,1,self.nobs]),axis=0) 

                # For BPDF
                BPDFparam = self.BPDFparam[:,:,self.iGood]

                # Loop through one by one
                # Some land covers do not have polarization (i.e. urban)
                I = np.zeros([self.nobs,1])
                Q = np.zeros([self.nobs,1])
                U = np.zeros([self.nobs,1])
                reflectance = np.zeros([self.nobs,1])
                surf_reflectance = np.zeros([self.nobs,1])
                ROT = np.zeros([nlev,self.nobs,1])

                for p in range(self.nobs):
                    if BPDFparam[2,0,p] == MISSING: 
                        args = [self.channel,
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
                        I_, reflectance_, ROT_, surf_reflectance_, Q_, U_, rc = BRDFvlidortWrapper(*args)                         

                    else:
                        args = [self.channel,
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
                        I_, reflectance_, ROT_, surf_reflectance_, Q_, U_, rc = vlidortWrapper(*args)
            
                    I[p:p+1,:] = I_
                    Q[p:p+1,:] = Q_
                    U[p:p+1,:] = U_
                    reflectance[p:p+1,:] = reflectance_
                    surf_reflectance[p:p+1,:] = surf_reflectance_
                    ROT[:,p:p+1,:] = ROT_
            
            elif self.albedoType == 'LAMBERTIAN':
                albedo = self.albedo[self.iGood,:]
                
                args = [self.channel,tau, ssa, pmom, 
                        pe, ze, te, 
                        albedo, 
                        sza, raa, vza, 
                        MISSING,
                        self.verbose]

                # Call VLIDORT wrapper function
                I, reflectance, ROT, Q, U, rc = vlidortWrapper(*args)  
                surf_reflectance = albedo
            if i == 0:
                self.ROT[self.iGood,:] = np.squeeze(ROT).T
            self.I[self.iGood,i] = np.squeeze(I)
            self.reflectance[self.iGood,i] = np.squeeze(reflectance)
            self.surf_reflectance[self.iGood,i] = np.squeeze(surf_reflectance)
            self.Q[self.iGood,i] = np.squeeze(Q)
            self.U[self.iGood,i] = np.squeeze(U) 


        self.writeNC()
        
    #---
    def writeNC (self,zlib=True):
        """
        Write a NetCDF file vlidort output
        """
        km = 72

        # Open NC file
        # ------------
        nc = Dataset(self.outFile,'w',format='NETCDF4_CLASSIC')

        # Set global attributes
        # ---------------------
        nc.title = 'VLIDORT Simulation of GEOS-5 multiangle polarized reflectance'
        nc.institution = 'NASA/Goddard Space Flight Center'
        nc.source = 'Global Model and Assimilation Office'
        nc.history = 'VLIDORT simulation run on sampled GEOS-5'
        nc.satangles  = self.VZA
        nc.references = 'n/a'
        nc.contact = 'Patricia Castellanos <patricia.castellanos@nasa.gov>'
        nc.Conventions = 'CF'
        nc.inFile = self.inFile
     
        # Create dimensions
        # -----------------
        nt = nc.createDimension('time',len(self.tyme))
        ls = nc.createDimension('ls',19)
        nz = nc.createDimension('lev',km)
        x  = nc.createDimension('x',1)
        y  = nc.createDimension('y',1)
        na = nc.createDimension('view_angles',2*len(self.VZA))

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


        vza = nc.createVariable('sensor_zenith','f4',('view_angles',),zlib=zlib)
        vza.long_name     = "sensor viewing zenith angle (VZA)"
        vza.missing_value = MISSING
        vza.units         = "degrees (positive forward view)"
        vza[:]            = np.append(-1.*self.VZA[::-1],self.VZA)

        vaa = nc.createVariable('sensor_azimuth','f4',('time','view_angles',),zlib=zlib)
        vaa.long_name     = "sensor viewing azimuth angle (VAA)"
        vaa.missing_value = MISSING
        vaa.units         = "degrees clockwise from North"
        
        vaaf = self.VAAf
        vaaf.shape = vaaf.shape + (1,)
        vaaf = np.repeat(vaaf,len(self.VZA),axis=1)
        
        
        vaab = self.VAAb
        vaab.shape = vaab.shape + (1,)        
        vaab = np.repeat(vaab,len(self.VZA),axis=1)

        vaa[:]            = np.append(vaab,vaaf,axis=1)

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
        ref = nc.createVariable('toa_reflectance','f4',('time','view_angles',),zlib=zlib,fill_value=MISSING)
        ref.standard_name = '%.2f nm TOA Reflectance' %self.channel
        ref.long_name     = '%.2f nm reflectance at the top of the atmosphere' %self.channel
        ref.missing_value = MISSING
        ref.units         = "None"
        ref[:]            = self.reflectance

        i = nc.createVariable('I','f4',('time','view_angles',),zlib=zlib,fill_value=MISSING)
        i.standard_name = '%.2f nm TOA I' %self.channel
        i.long_name     = '%.2f nm intensity at the top of the atmosphere' %self.channel
        i.missing_value = MISSING
        i.units         = "W m-2 sr-1 nm-1"
        i[:]            = self.I

        q = nc.createVariable('Q','f4',('time','view_angles',),zlib=zlib,fill_value=MISSING)
        q.standard_name = '%.2f nm TOA Q' %self.channel
        q.long_name     = '%.2f nm Q-component of the stokes vector at the top of the atmopshere' %self.channel
        q.missing_value = MISSING
        q.units         = "W m-2 sr-1 nm-1"
        q[:]            = self.Q    

        u = nc.createVariable('U','f4',('time','view_angles',),zlib=zlib,fill_value=MISSING)
        u.standard_name = '%.2f nm TOA U' %self.channel
        u.long_name     = '%.2f nm U-component of the stokes vector at the top of the atmopshere' %self.channel
        u.missing_value = MISSING
        u.units         = "W m-2 sr-1 nm-1"
        u[:]            = self.U

        sref = nc.createVariable('surf_reflectance','f4',('time','view_angles',),zlib=zlib,fill_value=MISSING)
        sref.standard_name = '%.2f nm Surface Reflectance' %self.channel
        sref.long_name     = '%.2f nm Bi-Directional Surface Reflectance' %self.channel
        sref.missing_value = MISSING
        sref.units         = "None"
        sref[:]            = self.surf_reflectance

        rot = nc.createVariable('ROT','f4',('time','lev',),zlib=zlib,fill_value=MISSING)
        rot.long_name = '%.2f nm Rayleigh Optical Thickness' %self.channel
        rot.missing_value = MISSING
        rot.units         = "None"
        rot[:]            = self.ROT

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

    inDir        = '/nobackup/3/pcastell/POLAR_LIDAR/CALIPSO/LevelB/Y{}/M{}'.format(date.year,str(date.month).zfill(2))
    inFile       = '{}/calipso-g5nr.lb2.%col.{}_{}z.nc4'.format(inDir,nymd,hour)
    albedoDir    = '/nobackup/3/pcastell/POLAR_LIDAR/CALIPSO/BRDF/MCD43C1/006/Y{}/M{}'.format(date.year,str(date.month).zfill(2))
    albedoFile   = '{}/calipso-g5nr.lb2.brdf.{}_{}z.nc4'.format(albedoDir,nymd,hour)
    ndviDir      = '/nobackup/3/pcastell/POLAR_LIDAR/CALIPSO/BPDF/NDVI/MYD13C2/006/Y{}/M{}'.format(date.year,str(date.month).zfill(2))
    ndviFile     = '{}/calipso-g5nr.lb2.ndvi.{}_{}z.nc4'.format(ndviDir,nymd,hour)
    lcDir        = '/nobackup/3/pcastell/POLAR_LIDAR/CALIPSO/BPDF/LAND_COVER/MCD12C1/051/Y{}/M{}'.format(date.year,str(date.month).zfill(2))
    lcFile       = '{}/calipso-g5nr.lb2.land_cover.{}_{}z.nc4'.format(lcDir,nymd,hour)    
    albedoType   = 'MODIS_BRDF_BPDF'

    channel  = 470
    chd      = get_chd(channel)
    outDir    = '/nobackup/3/pcastell/POLAR_LIDAR/CALIPSO/LevelC2/Y{}/M{}'.format(date.year,str(date.month).zfill(2))
    outFile   = '{}/calipso-g5nr.vlidort.vector.MCD43C.{}_{}z_{}nm.nc4'.format(outDir,nymd,hour,chd)
    
    rcFile   = 'Aod_EOS.rc'
    VZAname  = 'POLDER'
    orbit    = 'LEO'
    verbose  = True

    # Initialize VLIDORT class getting aerosol optical properties
    # -----------------------------------------------------------
    vlidort = POLAR_VLIDORT(inFile,outFile,rcFile,
                            albedoFile,albedoType,
                            channel,
                            VZAdic[VZAname],
                            HGTdic[orbit],
                            ndviFile=ndviFile,
                            lcFile=lcFile,
                            verbose=verbose)

   
    # Run ext_sampler
    # vlidort.runExt()

    # Run VLIDORT
    if vlidort.nobs > 0:
        vlidort.runVLIDORT()