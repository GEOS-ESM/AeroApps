#!/usr/bin/env python

"""
    Makes a polar plot of a typical scene - for benchmarking purposes

    Adapted from polar_vlidort
    Patricia Castellanos, May, 2017

"""

import os
from   netCDF4 import Dataset
from   mieobs  import getEdgeVars
import numpy   as np
from scipy import interpolate
from MAPL.constants import *
import VLIDORT_POLAR_ 


nMom     = 300

VZA = np.linspace(0,80,17)
VAA = np.linspace(0,360,73)
SZA = [30.0]
SAA = [0.0]
pixel = 1171

NDVI = 0.11846055
BPDFcoef = 5.99  

SurfaceFuncs = {'MODIS_BRDF'     : 'readSampledMODISBRDF',
                'MODIS_BRDF_BPDF': 'readSampledMODISBRDF',
                'LAMBERTIAN'     : 'readSampledLER'}

WrapperFuncs = {'MODIS_BRDF'     : VLIDORT_POLAR_.vector_brdf_modis,
                'MODIS_BRDF_BPDF': VLIDORT_POLAR_.vector_brdf_modis_bpdf,
                'BPDF'           : VLIDORT_POLAR_.vector_bpdf,
                'LAMBERTIAN'     : VLIDORT_POLAR_.vector_lambert}   

LandAlbedos  = 'MODIS_BRDF','MODIS_BRDF_BPDF','LAMBERTIAN','BPDF'

SDS_AER    = ['DELP','PS','RH','AIRDENS']
inDir      = '/discover/nobackup/projects/gmao/osse2/pub/c1440_NR/OBS/POLAR_LIDAR/CALIPSO/LevelB/Y2006/M08/'
inFile     = '{}/calipso-g5nr.lb2.aer_Nv.20060801_00z.nc4'.format(inDir)
inDir      = '/discover/nobackup/projects/gmao/osse2/pub/c1440_NR/OBS/POLAR_LIDAR/CALIPSO/LevelB/surface/BRDF/MCD43C1/006/Y2006/M08'
brdfFile   = '{}/calipso-g5nr.lb2.brdf.20060801_00z.nc4'.format(inDir)
inDir      = '/discover/nobackup/projects/gmao/osse2/pub/c1440_NR/OBS/POLAR_LIDAR/CALIPSO/LevelB/surface/SurfLER/Y2006/M08'
lerFile    = '{}/calipso-g5nr.lb2.omi_ler.20060801_00z.nc4'.format(inDir)
MISSING = -1.e+20

class POLAR_VLIDORT(object):
    """
    Everything needed for calling VLIDORT
    """
    def __init__(self,albedoType,channel,outFile,
                verbose=False):
        self.outFile = outFile
        self.albedoType = albedoType
        self.channel = channel
        self.verbose = verbose
        self.nMom    = nMom
        self.nobs    = 1

        self.SDS_AER = SDS_AER
        self.inFile  = inFile
        self.brdfFile = brdfFile
        self.lerFile  = lerFile

            
        # Read in surface data
        if (self.channel < 470) & ("MODIS_BRDF" in albedoType):
            albedoReader = getattr(self,'readHybridMODISBRDF')
        elif albedoType == 'BPDF':
            albedoReader = None
        else:
            albedoReader = getattr(self,SurfaceFuncs[albedoType])

        if albedoReader is not None:
            albedoReader()

        # Polarization
        if 'BPDF' in albedoType:
            self.BPDFinputs()

        # Calculate atmospheric profile properties needed for Rayleigh calc
        self.readSampledGEOS()
        self.computeAtmos()

        # Calculate Scene Geometry
        self.VZA = VZA
        self.VAA = VAA
        self.SZA = np.array(SZA)
        self.SAA = np.array(SAA)
        self.calcAngles()

    # --
    def BPDFinputs(self):
        """
        Set in NDVI and Landuse Coefficient
        """

        #BPDFparam(nparam,nch,nobs)
        # values for open shrubland
        self.BPDFparam = np.zeros([3,1,self.nobs])
        self.BPDFparam[0,0,:] = 1.5 
        self.BPDFparam[1,0,:] = NDVI
        self.BPDFparam[2,0,:] = BPDFcoef      


    # --
    def computeAtmos(self):

        pe, ze, te = getEdgeVars(self)

        self.pe = pe # (km,nobs)
        self.ze = ze 
        self.te = te

    # --
    def calcAngles(self):
        # define RAA according to photon travel direction
        saa = self.SAA + 180.0
        I = saa >= 360.
        saa[I] = saa[I] - 360.

        self.RAA = self.VAA - saa


    #---
    def readSampledGEOS(self):
        """
        Read in model sampled track
        """
        nc       = Dataset(self.inFile)

        for sds in self.SDS_AER:
            self.__dict__[sds] = []
            if len(nc.variables[sds][:].shape) == 1:
                var = np.array([nc.variables[sds][pixel]])
            else:
                var = nc.variables[sds][pixel,:]

            var.shape = (1,) + var.shape
            self.__dict__[sds] = var

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
            print 'opening BRDF file ',self.brdfFile
        nc = Dataset(self.brdfFile)

        for sds in SDS:
            self.__dict__[sds] = np.array([nc.variables[sds][pixel]])

        missing_value = nc.variables[sds].missing_value
        nc.close()

        nobs = self.nobs
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
            print 'opening BRDF abledo file ',self.brdfFile
        nc = Dataset(self.brdfFile)

        for sds in mSDS:
            self.__dict__[sds] = np.array([nc.variables[sds][pixel]])

        missing_value = nc.variables[sds].missing_value
        nc.close()

        if self.verbose:
            print 'opening LER albedo file ',self.lerFile
        nc = Dataset(self.lerFile)

        self.__dict__[lSDS] = np.array([nc.variables[lSDS][pixel]])
        nc.close()

        nobs = self.nobs
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


        self.albedo = np.array([0])  
        self.albedo.shape = (1,1)


    def runVLIDORT(self):
        """
        Calls VLIDORT 
        """
        pe   = self.pe
        ze   = self.ze
        te   = self.te

        nlev = pe.shape[0] - 1
        tau  = np.zeros([nlev,1,1])
        ssa  = tau
        pmom = np.zeros([nlev,1,1,nMom,6])
        

        nvza = len(self.VZA)
        nraa = len(self.RAA)
        sza  = self.SZA[0]

        self.I = np.ones([nvza,nraa])*MISSING
        self.Q = np.ones([nvza,nraa])*MISSING
        self.U = np.ones([nvza,nraa])*MISSING
        self.reflectance = np.ones([nvza,nraa])*MISSING
        self.surf_reflectance = np.ones([nvza,nraa])*MISSING
        self.BR_Q = np.ones([nvza,nraa])*MISSING
        self.BR_U = np.ones([nvza,nraa])*MISSING
        self.ROT = np.ones([nlev])*MISSING

        # Get VLIDORT wrapper name from dictionary
        vlidortWrapper = WrapperFuncs[self.albedoType]

        # loop through viewing angles and run VLIDORT
        for ivza,vza in enumerate(self.VZA):
            for iraa,raa in enumerate(self.RAA):
                if raa < 0:
                    raa = raa +360.0
                    
                # get args list for each surface model
                if self.albedoType == 'MODIS_BRDF':
                    kernel_wt = self.kernel_wt[:,:,0]
                    param     = self.RTLSparam[:,:,0]                
                    

                    args = [self.channel,tau, ssa, pmom, 
                            pe, ze, te, 
                            kernel_wt, param, 
                            sza, raa, vza, 
                            MISSING,
                            self.verbose]

                    # Call VLIDORT wrapper function
                    I, reflectance, ROT, surf_reflectance, Q, U, BR_Q, BR_U, rc = vlidortWrapper(*args)                        
                    
                elif self.albedoType == 'MODIS_BRDF_BPDF':
                    # For albedo
                    kernel_wt = self.kernel_wt[:,:,0]
                    RTLSparam = self.RTLSparam[:,:,0] 
                    RTLSparam = np.append(RTLSparam,np.zeros([1,1]),axis=0)

                    # For BPDF
                    BPDFparam = self.BPDFparam[:,:,0]

                    args = [self.channel,
                            tau, 
                            ssa, 
                            pmom, 
                            pe, 
                            ze, 
                            te, 
                            kernel_wt, 
                            RTLSparam, 
                            BPDFparam,
                            sza, 
                            raa, 
                            vza, 
                            MISSING,
                            self.verbose]

                    # Call VLIDORT wrapper function
                    I, reflectance, ROT, surf_reflectance, Q, U, BR_Q, BR_U, rc = vlidortWrapper(*args)

                elif self.albedoType == 'BPDF':

                    # For BPDF
                    BPDFparam = self.BPDFparam[:,:,0]

                    args = [self.channel,
                            tau, 
                            ssa, 
                            pmom, 
                            pe, 
                            ze, 
                            te, 
                            BPDFparam,
                            sza, 
                            raa, 
                            vza, 
                            MISSING,
                            self.verbose]

                    # Call VLIDORT wrapper function
                    I, reflectance, ROT, surf_reflectance, Q, U, BR_Q, BR_U, rc = vlidortWrapper(*args)
                                
                elif self.albedoType == 'LAMBERTIAN':
                    albedo = self.albedo
                    
                    args = [self.channel,tau, ssa, pmom, 
                            pe, ze, te, 
                            albedo, 
                            sza, raa, vza, 
                            MISSING,
                            self.verbose]

                    # Call VLIDORT wrapper function
                    I, reflectance, ROT, Q, U, rc = vlidortWrapper(*args)  
                    surf_reflectance = albedo
                if (ivza == 0) and (iraa == 0):
                    self.ROT = np.squeeze(ROT).T
                self.I[ivza,iraa] = np.squeeze(I)
                self.reflectance[ivza,iraa] = np.squeeze(reflectance)
                self.surf_reflectance[ivza,iraa] = np.squeeze(surf_reflectance)
                self.Q[ivza,iraa] = np.squeeze(Q)
                self.U[ivza,iraa] = np.squeeze(U) 
                if self.albedoType != 'LAMBERTIAN':
                    self.BR_Q[ivza,iraa] = np.squeeze(BR_Q)
                    self.BR_U[ivza,iraa] = np.squeeze(BR_U) 


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
        nc.title = 'VLIDORT Simulation of GEOS-5 multiangle polarized reflectance'
        nc.institution = 'NASA/Goddard Space Flight Center'
        nc.source = 'Global Model and Assimilation Office'
        nc.history = 'VLIDORT simulation run on sampled GEOS-5'
        nc.satangles  = self.VZA
        nc.references = 'n/a'
        nc.contact = 'Patricia Castellanos <patricia.castellanos@nasa.gov>'
        nc.Conventions = 'CF'
        nc.inFile = self.inFile
        nc.SZA    = self.SZA[0]
        nc.SAA    = self.SAA[0]

        if 'BPDF' in self.albedoType:
            nc.NDVI = NDVI
            nc.BPDFcoef = BPDFcoef

        if 'BRDF' in self.albedoType:
            nc.Riso = self.kernel_wt[0,0,0]
            nc.Rgeo = self.kernel_wt[1,0,0]
            nc.Rvol = self.kernel_wt[2,0,0]
     
        # Create dimensions
        # -----------------
        nt = nc.createDimension('vza',len(self.VZA))
        ls = nc.createDimension('vaa',len(self.VAA))
        nz = nc.createDimension('lev',km)


        vza = nc.createVariable('sensor_zenith','f4',('vza',),zlib=zlib)
        vza.long_name     = "sensor viewing zenith angle (VZA)"
        vza.missing_value = MISSING
        vza.units         = "degrees (positive forward view)"
        vza[:]            = self.VZA

        vaa = nc.createVariable('sensor_azimuth','f4',('vaa',),zlib=zlib)
        vaa.long_name     = "sensor viewing azimuth angle (VAA)"
        vaa.missing_value = MISSING
        vaa.units         = "degrees clockwise from North"
        
        vaa[:]            = self.VAA


        # Write VLIDORT Outputs
        # ---------------------
        ref = nc.createVariable('toa_reflectance','f4',('vza','vaa',),zlib=zlib,fill_value=MISSING)
        ref.standard_name = '%.2f nm TOA Reflectance' %self.channel
        ref.long_name     = '%.2f nm reflectance at the top of the atmosphere' %self.channel
        ref.missing_value = MISSING
        ref.units         = "None"
        ref[:]            = self.reflectance

        i = nc.createVariable('I','f4',('vza','vaa',),zlib=zlib,fill_value=MISSING)
        i.standard_name = '%.2f nm TOA I' %self.channel
        i.long_name     = '%.2f nm intensity at the top of the atmosphere' %self.channel
        i.missing_value = MISSING
        i.units         = "W m-2 sr-1 nm-1"
        i[:]            = self.I

        q = nc.createVariable('Q','f4',('vza','vaa',),zlib=zlib,fill_value=MISSING)
        q.standard_name = '%.2f nm TOA Q' %self.channel
        q.long_name     = '%.2f nm Q-component of the stokes vector at the top of the atmopshere' %self.channel
        q.missing_value = MISSING
        q.units         = "W m-2 sr-1 nm-1"
        q[:]            = self.Q    

        u = nc.createVariable('U','f4',('vza','vaa',),zlib=zlib,fill_value=MISSING)
        u.standard_name = '%.2f nm TOA U' %self.channel
        u.long_name     = '%.2f nm U-component of the stokes vector at the top of the atmopshere' %self.channel
        u.missing_value = MISSING
        u.units         = "W m-2 sr-1 nm-1"
        u[:]            = self.U

        sref = nc.createVariable('surf_reflectance','f4',('vza','vaa',),zlib=zlib,fill_value=MISSING)
        sref.standard_name = '%.2f nm Surface Reflectance' %self.channel
        sref.long_name     = '%.2f nm Bi-Directional Surface Reflectance' %self.channel
        sref.missing_value = MISSING
        sref.units         = "None"
        sref[:]            = self.surf_reflectance

        sref = nc.createVariable('surf_reflectance_Q','f4',('vza','vaa',),zlib=zlib,fill_value=MISSING)
        sref.standard_name = '%.2f nm Surface Reflectance Q' %self.channel
        sref.long_name     = '%.2f nm Bi-Directional Surface Reflectance Q' %self.channel
        sref.missing_value = MISSING
        sref.units         = "None"
        sref[:]            = self.BR_Q

        sref = nc.createVariable('surf_reflectance_U','f4',('vza','vaa',),zlib=zlib,fill_value=MISSING)
        sref.standard_name = '%.2f nm Surface Reflectance U' %self.channel
        sref.long_name     = '%.2f nm Bi-Directional Surface Reflectance U' %self.channel
        sref.missing_value = MISSING
        sref.units         = "None"
        sref[:]            = self.BR_U

        rot = nc.createVariable('ROT','f4',('lev',),zlib=zlib,fill_value=MISSING)
        rot.long_name = '%.2f nm Rayleigh Optical Thickness' %self.channel
        rot.missing_value = MISSING
        rot.units         = "None"
        rot[:]            = self.ROT

        # Close the file
        # --------------
        nc.close()

        if self.verbose:
            print " <> wrote %s"%(self.outFile)

    def runVLIDORTold(self):
        """
        Calls VLIDORT 
        """
        pe   = self.pe
        ze   = self.ze
        te   = self.te

        nlev = pe.shape[0] - 1
        tau  = np.zeros([nlev,1,1])
        ssa  = tau
        pmom = np.zeros([nlev,1,1,nMom,6])
        

        nvza = len(self.VZA)
        nraa = len(self.RAA)
        sza  = self.SZA[0]

        self.I = np.ones([nvza,nraa])*MISSING
        self.Q = np.ones([nvza,nraa])*MISSING
        self.U = np.ones([nvza,nraa])*MISSING
        self.reflectance = np.ones([nvza,nraa])*MISSING
        self.surf_reflectance = np.ones([nvza,nraa])*MISSING
        self.BR_Q = np.ones([nvza,nraa])*MISSING
        self.BR_U = np.ones([nvza,nraa])*MISSING
        self.ROT = np.ones([nlev])*MISSING

        # Get VLIDORT wrapper name from dictionary
        vlidortWrapper = WrapperFuncs[self.albedoType]

        # loop through viewing angles and run VLIDORT
        for ivza,vza in enumerate(self.VZA):
            for iraa,raa in enumerate(self.RAA):
                if raa < 0:
                    raa = raa +360.0
                    
                # get args list for each surface model
                if self.albedoType == 'MODIS_BRDF':
                    kernel_wt = self.kernel_wt[:,:,0]
                    param     = self.RTLSparam[:,:,0]                
                    

                    args = [self.channel,tau, ssa, pmom, 
                            pe, ze, te, 
                            kernel_wt, param, 
                            sza, raa, vza, 
                            MISSING,
                            self.verbose]

                    # Call VLIDORT wrapper function
                    I, reflectance, ROT, surf_reflectance, Q, U, BR_Q, BR_U, rc = vlidortWrapper(*args)                        
                    
                elif self.albedoType == 'MODIS_BRDF_BPDF':
                    # For albedo
                    kernel_wt = self.kernel_wt[:,:,0:1]
                    RTLSparam = self.RTLSparam[:,:,0:1] 
                    RTLSparam = np.append(RTLSparam,np.zeros([1,1,self.nobs]),axis=0) 

                    # For BPDF
                    BPDFparam = self.BPDFparam[:,:,0:1]

                    # Loop through one by one
                    # Some land covers do not have polarization (i.e. urban)
                    I = np.zeros([self.nobs,1])
                    Q = np.zeros([self.nobs,1])
                    U = np.zeros([self.nobs,1])
                    reflectance = np.zeros([self.nobs,1])
                    surf_reflectance = np.zeros([self.nobs,1])
                    BR_Q = np.zeros([self.nobs,1])
                    BR_U = np.zeros([self.nobs,1])
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
                                    sza, 
                                    raa, 
                                    vza, 
                                    MISSING,
                                    self.verbose]

                            BRDFvlidortWrapper = WrapperFuncs['MODIS_BRDF']
                            # Call VLIDORT wrapper function
                            I_, reflectance_, ROT_, surf_reflectance_, Q_, U_, BR_Q_, BR_U_, rc = BRDFvlidortWrapper(*args)                         

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
                                    sza, 
                                    raa, 
                                    vza, 
                                    MISSING,
                                    self.verbose]

                            # Call VLIDORT wrapper function
                            I_, reflectance_, ROT_, surf_reflectance_, Q_, U_, BR_Q_, BR_U_, rc = vlidortWrapper(*args)
                
                        I[p:p+1,:] = I_
                        Q[p:p+1,:] = Q_
                        U[p:p+1,:] = U_
                        reflectance[p:p+1,:] = reflectance_
                        surf_reflectance[p:p+1,:] = surf_reflectance_
                        BR_Q[p:p+1,:] = BR_Q_
                        BR_U[p:p+1,:] = BR_U_
                        ROT[:,p:p+1,:] = ROT_
                
                elif self.albedoType == 'LAMBERTIAN':
                    albedo = self.albedo
                    
                    args = [self.channel,tau, ssa, pmom, 
                            pe, ze, te, 
                            albedo, 
                            sza, raa, vza, 
                            MISSING,
                            self.verbose]

                    # Call VLIDORT wrapper function
                    I, reflectance, ROT, Q, U, rc = vlidortWrapper(*args)  
                    surf_reflectance = albedo
                if (ivza == 0) and (iraa == 0):
                    self.ROT = np.squeeze(ROT).T
                self.I[ivza,iraa] = np.squeeze(I)
                self.reflectance[ivza,iraa] = np.squeeze(reflectance)
                self.surf_reflectance[ivza,iraa] = np.squeeze(surf_reflectance)
                self.Q[ivza,iraa] = np.squeeze(Q)
                self.U[ivza,iraa] = np.squeeze(U) 
                if self.albedoType != 'LAMBERTIAN':
                    self.BR_Q[ivza,iraa] = np.squeeze(BR_Q)
                    self.BR_U[ivza,iraa] = np.squeeze(BR_U) 


        self.writeNC()
    

def get_chd(channel):
    chd = '%.2f'%channel
    chd = chd.replace('.','d')

    return chd
    
#------------------------------------ M A I N ------------------------------------

if __name__ == "__main__":
    format   = 'NETCDF4_CLASSIC'

    rootDir  = '/discover/nobackup/projects/gmao/osse2/pub/c1440_NR/OBS/POLAR_LIDAR/CALIPSO/'
    albedoType   = 'BRDF_BPDF'

    channel  = 865
    chd      = get_chd(channel)
    outDir    = './benchmark_mueller_fix'
    outFile   = '{}/calipso-g5nr.vlidort.vector.BPDF.{}.nc4'.format(outDir,chd)
    
    rcFile   = 'Aod_EOS.rc'
    verbose  = True

    # Initialize VLIDORT class getting aerosol optical properties
    # -----------------------------------------------------------
    vlidort = POLAR_VLIDORT(albedoType,channel,outFile,
                            verbose=verbose)

   
    vlidort.runVLIDORT()
