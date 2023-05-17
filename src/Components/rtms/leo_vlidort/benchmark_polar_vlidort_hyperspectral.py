#!/usr/bin/env python3

"""
    Makes a polar plot of a typical scene - for benchmarking purposes

    Adapted from polar_vlidort
    Patricia Castellanos, May, 2017

"""

import os
from   netCDF4 import Dataset
from   mieobs  import getAOPvector, getEdgeVars, getAOPint, getAOPext
import numpy   as np
from scipy import interpolate
from MAPL.constants import *
import VLIDORT_POLAR_ 
from polar_vlidort import POLAR_VLIDORT
from scipy.special.orthogonal import legendre


nMom     = 300
nPol     = 6
km       = 72

VZA = np.linspace(0,80,17)
#VZA = np.linspace(0,80,9)
VAA = np.linspace(0,360,73)
#VAA = np.linspace(0,360,37)
SZA = [30.0]
SAA = [0.0]
pixel = 1171

NDVI = 0.11846055
BPDFcoef = 5.99


SurfaceFuncs = {'MODIS_BRDF'     : 'readSampledMODISBRDF',
                'MODIS_BRDF_BPDF': 'readSampledMODISBRDF',
                'LAMBERTIAN'     : 'readSampledLER',
                'LAMBERTIAN_BPDF': 'readSampledLER',
                'GissCX'         : 'readSampledWindCX',
                'OCIGissCX'      : 'readSampledWindCX',
                'OCIGissCX_NOBM_CLOUD'      : 'readSampledWindCX',
                'CX'             : 'readSampledWindCX',
                'OCICX'          : 'readSampledWindCX'}

WrapperFuncs = {'MODIS_BRDF'                : VLIDORT_POLAR_.vector_brdf_modis,
                'MODIS_BRDF_BPDF'           : VLIDORT_POLAR_.vector_brdf_modis_bpdf,
                'BPDF'                      : VLIDORT_POLAR_.vector_bpdf,
                'LAMBERTIAN'                : VLIDORT_POLAR_.vector_lambert,
                'LAMBERTIAN_BPDF'           : VLIDORT_POLAR_.vector_lambert_bpdf,
                'GissCX'                    : VLIDORT_POLAR_.vector_gisscx,
                'CX'                        : VLIDORT_POLAR_.vector_cx,
                'OCICX'                     : VLIDORT_POLAR_.vector_ocicx,
                'OCIGissCX'                 : VLIDORT_POLAR_.vector_ocigisscx,
                'OCIGissCX_NOBM_CLOUD'      : VLIDORT_POLAR_.vector_ocigisscx_nobm_cloud,
                'ROT_CALC'                  : VLIDORT_POLAR_.rot_calc}   


# Generic Lists of Varnames and Units
VNAMES_DU = ['DU001','DU002','DU003','DU004','DU005']
VNAMES_SS = ['SS001','SS002','SS003','SS004','SS005']
VNAMES_BC = ['BCPHOBIC','BCPHILIC']
VNAMES_OC = ['OCPHOBIC','OCPHILIC']
VNAMES_SU = ['SO4']

META    = ['DELP','PS','RH','AIRDENS','trjLon','trjLat','isotime']
AERNAMES = VNAMES_SU #+ VNAMES_SS + VNAMES_OC + VNAMES_BC + VNAMES_DU
#AERNAMES = ['DU001']
SDS_AER = META + AERNAMES
MISSING = -1.e+20

class BENCHMARK(POLAR_VLIDORT):
    """
    Everything needed for calling VLIDORT
    """
    def __init__(self,albedoType,channels,outFile,rootDir,rcFile,
                verbose=False,aerosol=True,dark=False,simple_aerosol=False,
                emde=False,grasp=False,grasp_higher=False,rayleigh=False,
                nstreams=12,plane_parallel=True,
                sleave_adjust=False):

        self.outFile = outFile
        self.albedoType = albedoType
        self.channels = channels
        self.verbose = verbose
        self.nMom    = nMom
        self.nPol    = nPol
        self.nobs    = 1
        self.nch     = len(channels)
        self.km      = km
        self.SDS_AER = SDS_AER
        self.SDS_META= META
        self.AERNAMES= AERNAMES
        self.rootDir = rootDir
        self.rcFile  = rcFile
        self.aerosol = aerosol
        self.simple_aerosol = simple_aerosol
        self.dark    = dark
        self.emde    = emde
        self.grasp   = grasp
        self.grasp_higher = grasp_higher
        self.rayleigh = rayleigh
        self.nstreams = nstreams
        self.plane_parallel = plane_parallel
        self.sleave_adjust = sleave_adjust

        # get filesnames
        self.getFilenames()
            
        # Read in surface data
        self.getSurface()

        # Polarization
        if 'BPDF' in albedoType:
            self.BPDFinputs()

        # Water leaving radiance
        if 'NOBM' in albedoType:
            self.readRRS()

        # Read Model Data
        self.readSampledGEOS()
        # Calculate atmospheric profile properties needed for Rayleigh calc      
        self.computeAtmos()
        # Calculate aerosol optical properties
        # set to zero if aerosol is false
        self.computeMie()          

        # Calculate Scene Geometry
        self.VZA = VZA
        self.VAA = VAA
        self.SZA = np.array(SZA)
        self.SAA = np.array(SAA)
        self.calcAngles()

    #---
    def getSurface(self):
        if 'CX' in self.albedoType:
            albedoReader = getattr(self,SurfaceFuncs[albedoType])
            albedoReader()
        elif albedoType == 'BPDF':
            pass
        else: 
            for ch in self.channels:
                if (ch < 470) & ("MODIS_BRDF" in self.albedoType):
                    albedoReader = getattr(self,'readHybridMODISBRDF')
                else:
                    albedoReader = getattr(self,SurfaceFuncs[albedoType])
                albedoReader()

    #---
    def getFilenames(self):
        inDir           = '{}/LevelB/Y2006/M08/'.format(self.rootDir)
        self.inFile     = '{}/calipso-g5nr.lb2.aer_Nv.20060801_00z.nc4'.format(inDir)
        inDir           = '{}/LevelB/surface/BRDF/MCD43C1/006/Y2006/M08'.format(self.rootDir)
        self.brdfFile   = '{}/calipso-g5nr.lb2.brdf.20060801_00z.nc4'.format(inDir)
        inDir           = '{}/LevelB/surface/SurfLER/Y2006/M08'.format(self.rootDir)
        self.lerFile    = '{}/calipso-g5nr.lb2.omi_ler.20060801_00z.nc4'.format(inDir)
        self.rrsFile    = '/nobackup/PACE/LevelB/surface/SLEAVE/NOBM/Y2006/M03/D24/pace-g5nr.lb.sleave.20060324_005000.nc4'

    #---
    def computeMie(self):
        """
        Computes aerosol optical quantities 
        """
        self.tau = []
        self.ssa = []
        self.g   = []
        self.pmom = []
        for ch in self.channels:
            rcFile = self.rcFile.format(ch)
            if self.simple_aerosol and self.aerosol:
                self.RH[:] = 0.00

                tauO,ssa,g,pmom = getAOPvector(self,[ch],
                                 vnames=self.AERNAMES,
                                 Verbose=True,
                                 rcfile=rcFile,
                                 nMom=self.nMom)

                dz = self.ze[:-1,0] - self.ze[1:,0]  
                z  = self.ze[:-1,0] - 0.5*dz
                H  = 1000 # m
                ext = np.exp(-z/H)
                # scale so integral equals AOD
                total = np.sum(ext*dz)
                ext   = ext*tauO.sum()/total
                tau   = ext*dz  
                tau.shape = tau.shape + (1,1)
                tau       = tau

                # i think this sign needs to be changed
                pmom[:,:,:,:,1] = -1.*pmom[:,:,:,:,1]
                pmom[:,:,:,:,3] = -1.*pmom[:,:,:,:,3]

            elif self.aerosol:
                tau,ssa,g,pmom = getAOPvector(self,[ch],
                                 vnames=self.AERNAMES,
                                 Verbose=True,
                                 rcfile=rcFile,
                                 nMom=self.nMom)
            else:
                tau = np.zeros([self.km,self.nch,self.nobs])
                ssa = np.zeros([self.km,self.nch,self.nobs])
                g   = np.zeros([self.km,self.nch,self.nobs])
                pmom = np.zeros([self.km,self.nch,self.nobs,self.nMom,self.nPol])

            self.tau.append(tau)  #(km,nch,nobs)
            self.ssa.append(ssa)  #(km,nch,nobs)
            self.g.append(g)    #(km,nch,nobs)
            self.pmom.append(pmom)  #(km,nch,nobs,nMom,nPol)
        
        self.tau  = np.concatenate(self.tau,axis=1)
        self.ssa  = np.concatenate(self.ssa,axis=1)
        self.g    = np.concatenate(self.g,axis=1)
        self.pmom = np.concatenate(self.pmom,axis=1)

    # --
    def readSampledWindCX(self):
        """
        Cox-munk inputs
        """
        self.mr = np.repeat([1.333],self.nch)
        self.U10m = np.array([3.0])
        self.V10m = np.array([4.0])

    def readRRS(self):
        """
        sun normalized water leaving radiance
        divided by extraterrestrial solar irradiance
        this has to be multiplied by cos(SZA) to get
        the VLIDORT SL_ISOTROPIC value
        valid for 400-450 nm
        """

        nc = Dataset(self.rrsFile)
        rrs = nc.variables['rrs'][0,:,900,600]
        wav = nc.variables['wavelength'][:]
        f = interpolate.interp1d(wav,rrs)
        self.RRS = f(self.channels)
        i = self.RRS<0
        self.RRS[i] = 0.0

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
        if self.aerosol:
            SDS = self.SDS_AER
        else:
            SDS = self.SDS_META
        for sds in SDS:
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
            chmin = np.argmin(dch[dch>0])
            chmin = MODIS_channels[dch>0][chmin]
            chmax = np.argmax(dch[dch<0])
            chmax = MODIS_channels[dch<0][chmax]

            SDS = 'Riso'+chmin,'Rgeo'+chmin,'Rvol'+chmin,'Riso'+chmax,'Rgeo'+chmax,'Rvol'+chmax

        if self.verbose:
            print('opening BRDF file ',self.brdfFile)
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
            print('opening BRDF abledo file ',self.brdfFile)
        nc = Dataset(self.brdfFile)

        for sds in mSDS:
            self.__dict__[sds] = np.array([nc.variables[sds][pixel]])

        missing_value = nc.variables[sds].missing_value
        nc.close()

        if self.verbose:
            print('opening LER albedo file ',self.lerFile)
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

    #---
    def readSampledLER(self):

        if self.dark:
            self.albedo = np.array([0.0])  
        else:
            self.albedo = np.array([0.3])  

        self.albedo.shape = (1,1)

    #---
    def rot2scatplane(self):
        """
        Rotate meridonal plane (VLIDORT) outputs to scattering plane (GRASP)
        adapted from Kirk's code rot2scatplane.pro
        """
       
        nch  = self.nch
        nvza = len(self.VZA)
        nraa = len(self.VAA)

        #raa = np.radians(self.VAA - self.SAA)
        raa = np.radians(self.RAA)
        razi_pm = np.ones(nraa)
        razi_pm[np.sin(raa) >= 0] = -1.0
        cos_razi = np.cos(raa)

        u_s = np.cos(np.radians(self.SZA[0]))
        u_v = np.cos(np.radians(self.VZA))

        root_s = np.sqrt(1.0-u_s**2)
        self.Q_out = np.zeros([nch,nvza,nraa])
        self.U_out = np.zeros([nch,nvza,nraa])
        self.Qsurf_out = np.zeros([nch,nvza,nraa])
        self.Usurf_out = np.zeros([nch,nvza,nraa])        
        for ivza in range(nvza):
            for iraa in range(nraa):
                #print 'ivza,iraa',ivza,iraa
                root_v = np.sqrt(1.0-u_v[ivza]**2)
         
                cos_theta= -1.0*u_s*u_v[ivza] + root_v*root_s*cos_razi[iraa]
                root_theta= np.sqrt(1.0 - cos_theta**2)
                scatang = np.arccos(cos_theta) #scattering angle for output (not used below)
        
                # equation 3.16 in Hansen and Travis --------------------------
                # Special limit case --------------------
                #if np.abs(cos_theta) > 0.999999:
                if np.abs(cos_theta) == 1.0:
                    cosi1 = 0.0
                    cosi2 = 0.0
                else:
                    cosi1= razi_pm[iraa]*(-1.0*u_v[ivza]*root_s - u_s*root_v*cos_razi[iraa])/root_theta
                    cosi2= razi_pm[iraa]*(u_s*root_v + u_v[ivza]*root_s*cos_razi[iraa])/root_theta

                # equation (10) in Hovenier ------------------------
                # error correction for high cos(i)^2 values ------------
                #if (cosi1**2.0) > 0.999999:
                if (cosi1**2.0) >= 1.0:
                    cos2i1 = 1.0
                    sin2i1 = 0.0
                else:
                    sin2i1= 2.0*np.sqrt(1.0-(cosi1**2.0))*cosi1
                    cos2i1= 2.0*(cosi1**2.0)-1.0

                #if (cosi2**2.0) > 0.999999:
                if (cosi2**2.0) >= 1.0:
                    cos2i2 = 1.0
                    sin2i2 = 0.0
                else:
                    sin2i2= 2.0*np.sqrt(1.0-(cosi2**2.0))*cosi2
                    cos2i2= 2.0*(cosi2**2.0)-1.0

                # rotate into scattering plane as shown in (2) of Hovenier 
                q_in = self.Q[:,ivza,iraa]
                u_in = self.U[:,ivza,iraa]
                q_out= q_in*cos2i2 - u_in*sin2i2
                u_out= q_in*sin2i2 + u_in*cos2i2

                self.Q_out[:,ivza,iraa] = -1.*q_out
                self.U_out[:,ivza,iraa] = -1.*u_out

                q_in = self.BR_Q[:,ivza,iraa]
                u_in = self.BR_U[:,ivza,iraa]
                q_out= q_in*cos2i2 - u_in*sin2i2
                u_out= q_in*sin2i2 + u_in*cos2i2                
                self.Qsurf_out[:,ivza,iraa] = -1.*q_out
                self.Usurf_out[:,ivza,iraa] = -1.*u_out

    def A(self,X):
        """
        Rayleigh scattering cross section (cm2/molecule)
        """
        pi = np.pi
        Av  = 6.0221367e23    #Avogadros number
        Ts = 288.15           #K, standard temperature

        XX=1./X
        XX2=XX*XX
        CM = X/10000.0  # lambda in cm
        pi3 = pi*pi*pi

        #  (ns-1)*1e8  - refractive index of dry air at standard pressure and temperature for CO2 = 300 ppmv
        RI = 8060.51 + 2480990.0/(132.274-XX2) + 17455.7/(39.32957-XX2)
        RI = RI*1e-8 + 1.0

        C = 360.0*1e-6 # parts per volume
        RI = (RI-1)*(1.0 + 0.54*(C - 0.0003)) + 1
        RI2 = RI*RI

        # adjust for CO2 = 360 ppm
        Ns = (Av/22.4141)*(273.15/Ts)*(1.0/1000)
        A = 24.0*pi3*self.F_air(XX2,C)* ((RI2-1)*(RI2-1))/(Ns*Ns*CM*CM*CM*CM*(RI2+2)*(RI2+2))

        return A

    def F_air(self,XX2,C):
        """
        King factor - depolarization
        """
        depol1 = 1.034 + 3.17E-4 * XX2  # N2
        depol2 = 1.096 + 1.385E-3 * XX2 + 1.448E-4 * XX2 * XX2  # O2
        depol3 = 1.  # Ar
        depol4 = 1.15  # CO2

        # C*100  conc in percent parts per volume
        F_air= (78.084*depol1 + 20.946*depol2 + 0.934*depol3 +  C*100.0*depol4) / (78.084 + 20.946 + 0.934 + C*100.0)

        return F_air


    def runVLIDORT(self):
        """
        Calls VLIDORT 
        """
        pe   = self.pe
        ze   = self.ze
        te   = self.te

        nlev = self.km
        tau  = self.tau   
        ssa  = self.ssa   
        pmom = self.pmom  
        
        tauC = self.tau*0.0
        ssaC = self.ssa*0.0
        pmomC = self.pmom*0.0

        nch = self.nch

        nvza = len(self.VZA)
        nraa = len(self.RAA)
        sza  = self.SZA[0]

        sleave = self.RRS*np.cos(np.radians(sza))

        self.I = np.ones([nch,nvza,nraa])*MISSING
        self.Q = np.ones([nch,nvza,nraa])*MISSING
        self.U = np.ones([nch,nvza,nraa])*MISSING
        self.reflectance = np.ones([nch,nvza,nraa])*MISSING
        self.surf_reflectance = np.ones([nch,nvza,nraa])*MISSING
        self.BR_Q = np.ones([nch,nvza,nraa])*MISSING
        self.BR_U = np.ones([nch,nvza,nraa])*MISSING
        self.adjusted_sleave = np.ones([nch,nvza,nraa])*MISSING

        # Calculate ROT
        args = [self.channels, pe, ze, te, MISSING, self.verbose]
        vlidortWrapper = WrapperFuncs['ROT_CALC']
        ROT, depol_ratio, rc = vlidortWrapper(*args) 
        
        self.ROT = np.squeeze(ROT).T
        alpha = ROT*0.0

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
                    

                    args = [self.channels, self.nstreams, self.plane_parallel, ROT, depol_ratio, tau, ssa, pmom, 
                            pe, ze, te, 
                            kernel_wt, param, 
                            sza, raa, vza, 
                            MISSING,
                            self.verbose]

                    # Call VLIDORT wrapper function
                    I, reflectance, surf_reflectance, Q, U, BR_Q, BR_U, rc = vlidortWrapper(*args)                        
                    
                elif self.albedoType == 'MODIS_BRDF_BPDF':
                    # For albedo
                    kernel_wt = self.kernel_wt[:,:,0]
                    RTLSparam = self.RTLSparam[:,:,0] 
                    RTLSparam = np.append(RTLSparam,np.zeros([1,1]),axis=0)

                    # For BPDF
                    BPDFparam = self.BPDFparam[:,:,0]

                    args = [self.channels,
                            self.nstreams,
                            self.plane_parallel, 
                            ROT,
                            depol_ratio, 
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
                    I, reflectance, surf_reflectance, Q, U, BR_Q, BR_U, rc = vlidortWrapper(*args)

                elif self.albedoType == 'BPDF':

                    # For BPDF
                    BPDFparam = self.BPDFparam[:,:,0]

                    args = [self.channels,
                            self.nstreams,
                            self.plane_parallel, 
                            ROT,
                            depol_ratio, 
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
                    I, reflectance, surf_reflectance, Q, U, BR_Q, BR_U, rc = vlidortWrapper(*args)

                elif self.albedoType == 'LAMBERTIAN_BPDF':
                    albedo = self.albedo
                    # For BPDF
                    BPDFparam = self.BPDFparam[:,:,0]

                    args = [self.channels,
                            self.nstreams,
                            self.plane_parallel, 
                            ROT,
                            depol_ratio, 
                            tau, 
                            ssa, 
                            pmom, 
                            pe, 
                            ze, 
                            te, 
                            albedo,
                            BPDFparam,
                            sza, 
                            raa, 
                            vza, 
                            MISSING,
                            self.verbose]

                    # Call VLIDORT wrapper function
                    I, reflectance, surf_reflectance, Q, U, BR_Q, BR_U, rc = vlidortWrapper(*args)                    
                
                elif ('CX' in self.albedoType) and (self.albedoType != 'OCIGissCX_NOBM_CLOUD'):
                    args = [self.channels, self.nstreams, self.plane_parallel, ROT, depol_ratio, tau, ssa, pmom,
                            pe, ze, te,
                            self.U10m, self.V10m, self.mr,
                            sza, raa, vza,
                            MISSING,
                            self.verbose]

                    I, reflectance, surf_reflectance, Q, U, BR_Q, BR_U, rc = vlidortWrapper(*args)
                elif self.albedoType == 'OCIGissCX_NOBM_CLOUD':
                    args = [self.channels, self.nstreams, 
                            self.plane_parallel, ROT, depol_ratio, alpha,
                            tau, ssa, pmom,
                            tauC, ssaC, pmomC,
                            tauC, ssaC, pmomC,
                            pe, ze, te,
                            self.U10m, self.V10m, self.mr,
                            sleave, self.sleave_adjust,
                            sza, raa, vza,
                            MISSING,
                            self.verbose]
                    
                    I, reflectance, surf_reflectance, Q, U, BR_Q, BR_U, rc, adjusted_sleave = vlidortWrapper(*args)                    
                elif self.albedoType == 'LAMBERTIAN':
                    albedo = self.albedo
                    
                    args = [self.channel, self.nstreams, self.plane_parallel, ROT, depol_ratio, alpha, tau, ssa, pmom, 
                            pe, ze, te, 
                            albedo, 
                            sza, raa, vza, 
                            MISSING,
                            self.verbose]

                    # Call VLIDORT wrapper function
                    I, reflectance, Q, U, rc = vlidortWrapper(*args)  
                    surf_reflectance = albedo

                self.I[:,ivza,iraa] = np.squeeze(I)
                self.reflectance[:,ivza,iraa] = np.squeeze(reflectance)
                self.surf_reflectance[:,ivza,iraa] = np.squeeze(surf_reflectance)
                self.Q[:,ivza,iraa] = np.squeeze(Q)
                self.U[:,ivza,iraa] = np.squeeze(U) 
                if self.albedoType != 'LAMBERTIAN':
                    self.BR_Q[:,ivza,iraa] = np.squeeze(BR_Q)
                    self.BR_U[:,ivza,iraa] = np.squeeze(BR_U)
                if ('NOBM' in self.albedoType) and (self.sleave_adjust is True):
                    self.adjusted_sleave[:,ivza,iraa] = np.squeeze(adjusted_sleave)


        self.rot2scatplane()
        self.writeNC()
        
    #---
    def writeNC (self,zlib=True):
        """
        Write a NetCDF file vlidort output
        """

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
        ch = nc.createDimension('nch',self.nch)
        nt = nc.createDimension('vza',len(self.VZA))
        ls = nc.createDimension('vaa',len(self.VAA))
        nz = nc.createDimension('lev',km)
        nze = nc.createDimension('leve',km+1)
        no = nc.createDimension('nobs',1)   
        # if self.aerosol:
        #     nr = nc.createDimension('radius',len(self.R))  
        #     na = nc.createDimension('angle',len(self.angle))    
        ch = nc.createVariable('channels','f4',('nch',),zlib=zlib)
        ch.long_name = "wavelength"
        ch.units     = "nm"
        ch[:]        = self.channels

        vza = nc.createVariable('sensor_zenith','f4',('vza',),zlib=zlib)
        vza.long_name     = "sensor viewing zenith angle (VZA)"
        vza.missing_value = np.float32(MISSING)
        vza.units         = "degrees (positive forward view)"
        vza[:]            = self.VZA

        vaa = nc.createVariable('sensor_azimuth','f4',('vaa',),zlib=zlib)
        vaa.long_name     = "sensor viewing azimuth angle (VAA)"
        vaa.missing_value = np.float32(MISSING)
        vaa.units         = "degrees clockwise from North"
        
        vaa[:]            = self.VAA


        # Write VLIDORT Outputs
        # ---------------------
        ref = nc.createVariable('toa_reflectance','f4',('nch','vza','vaa',),zlib=zlib,fill_value=MISSING)
        ref.standard_name = 'spectral TOA Reflectance'
        ref.long_name     = 'spectral reflectance at the top of the atmosphere'
        ref.missing_value = np.float32(MISSING)
        ref.units         = "None"
        ref[:]            = self.reflectance

        i = nc.createVariable('I','f4',('nch','vza','vaa',),zlib=zlib,fill_value=MISSING)
        i.standard_name = 'spectral TOA I'
        i.long_name     = 'spectral intensity at the top of the atmosphere'
        i.missing_value = np.float32(MISSING)
        i.units         = "W m-2 sr-1 nm-1"
        i[:]            = self.I

        q = nc.createVariable('Q','f4',('nch','vza','vaa',),zlib=zlib,fill_value=MISSING)
        q.standard_name = 'spectral TOA Q'
        q.long_name     = 'spectral Q-component of the stokes vector at the top of the atmopshere' 
        q.missing_value = np.float32(MISSING)
        q.units         = "W m-2 sr-1 nm-1"
        q[:]            = self.Q    

        u = nc.createVariable('U','f4',('nch','vza','vaa',),zlib=zlib,fill_value=MISSING)
        u.standard_name = 'spectral TOA U'
        u.long_name     = 'spectral U-component of the stokes vector at the top of the atmopshere' 
        u.missing_value = np.float32(MISSING)
        u.units         = "W m-2 sr-1 nm-1"
        u[:]            = self.U

        q = nc.createVariable('Q_scatplane','f4',('nch','vza','vaa',),zlib=zlib,fill_value=MISSING)
        q.standard_name = 'spectral TOA Q rotated to scattering plane'
        q.long_name     = 'spectral Q-component of the stokes vector at the top of the atmopshere'
        q.missing_value = np.float32(MISSING)
        q.units         = "W m-2 sr-1 nm-1"
        q[:]            = self.Q_out

        u = nc.createVariable('U_scatplane','f4',('nch','vza','vaa',),zlib=zlib,fill_value=MISSING)
        u.standard_name = 'spectral TOA U rotated to scattering plane'
        u.long_name     = 'spectral U-component of the stokes vector at the top of the atmopshere' 
        u.missing_value = np.float32(MISSING)
        u.units         = "W m-2 sr-1 nm-1"
        u[:]            = self.U_out
        sref = nc.createVariable('surf_reflectance','f4',('nch','vza','vaa',),zlib=zlib,fill_value=MISSING)
        sref.standard_name = 'spectral Surface Reflectance' 
        sref.long_name     = 'spectral Bi-Directional Surface Reflectance'
        sref.missing_value = np.float32(MISSING)
        sref.units         = "None"
        sref[:]            = self.surf_reflectance

        sref = nc.createVariable('surf_reflectance_Q','f4',('nch','vza','vaa',),zlib=zlib,fill_value=MISSING)
        sref.standard_name = 'spectral Surface Reflectance Q'
        sref.long_name     = 'spectral Bi-Directional Surface Reflectance Q'
        sref.missing_value = np.float32(MISSING)
        sref.units         = "None"
        sref[:]            = self.BR_Q

        sref = nc.createVariable('surf_reflectance_U','f4',('nch','vza','vaa',),zlib=zlib,fill_value=MISSING)
        sref.standard_name = 'spectral Surface Reflectance U'
        sref.long_name     = 'spectral Bi-Directional Surface Reflectance U'
        sref.missing_value = np.float32(MISSING)
        sref.units         = "None"
        sref[:]            = self.BR_U

        sref = nc.createVariable('surf_reflectance_Q_scatplane','f4',('nch','vza','vaa',),zlib=zlib,fill_value=MISSING)
        sref.standard_name = 'spectral Surface Reflectance Q'
        sref.long_name     = 'spectral Bi-Directional Surface Reflectance Q'
        sref.missing_value = np.float32(MISSING)
        sref.units         = "None"
        sref[:]            = self.Qsurf_out

        sref = nc.createVariable('surf_reflectance_U_scatplane','f4',('nch','vza','vaa',),zlib=zlib,fill_value=MISSING)
        sref.standard_name = 'spectral Surface Reflectance U'
        sref.long_name     = 'spectral Bi-Directional Surface Reflectance U'
        sref.missing_value = np.float32(MISSING)
        sref.units         = "None"
        sref[:]            = self.Usurf_out


        rot = nc.createVariable('ROT','f4',('lev','nch',),zlib=zlib,fill_value=MISSING)
        rot.long_name = 'spectral Rayleigh Optical Thickness'
        rot.missing_value = np.float32(MISSING)
        rot.units         = "None"
        rot[:]            = self.ROT

        pe = nc.createVariable('PE','f4',('leve',),zlib=zlib,fill_value=MISSING)
        pe.long_name = 'layer edge pressure' 
        pe.missing_value = np.float32(MISSING)
        pe.units         = "Pa"
        pe[:]            = self.pe

        ze = nc.createVariable('ZE','f4',('leve',),zlib=zlib,fill_value=MISSING)
        ze.long_name = 'layer edge height' 
        ze.missing_value = np.float32(MISSING)
        ze.units         = "m"
        ze[:]            = self.ze

        te = nc.createVariable('TE','f4',('leve',),zlib=zlib,fill_value=MISSING)
        te.long_name = 'layer edge temperature' 
        te.missing_value = np.float32(MISSING)
        te.units         = "K"
        te[:]            = self.te

        dens = nc.createVariable('AIRDENS','f4',('lev',),zlib=zlib,fill_value=MISSING)
        dens.long_name = 'layer air density' 
        dens.missing_value = np.float32(MISSING)
        dens.units         = "kg m-3"
        dens[:]            = self.AIRDENS[0,:]      

        if 'BRDF' in self.albedoType:
            chs = str(int(self.channel))
            riso = nc.createVariable('RISO','f4',('nobs',),zlib=zlib,fill_value=MISSING)
            riso.long_name = 'RTLS isotropic kernel weight' 
            riso.missing_value = np.float32(MISSING)
            riso.units         = "None"
            riso[:]            = self.__dict__['Riso'+chs]

            rgeo = nc.createVariable('RGEO','f4',('nobs',),zlib=zlib,fill_value=MISSING)
            rgeo.long_name = 'RTLS geometric kernel weight' 
            rgeo.missing_value = np.float32(MISSING)
            rgeo.units         = "None"
            rgeo[:]            = self.__dict__['Rgeo'+chs]

            rvol = nc.createVariable('RVOL','f4',('nobs',),zlib=zlib,fill_value=MISSING)
            rvol.long_name = 'RTLS volumetric kernel weight' 
            rvol.missing_value = np.float32(MISSING)
            rvol.units         = "None"
            rvol[:]            = self.__dict__['Rvol'+chs]

        if 'CX' in self.albedoType:
            mr = nc.createVariable('mr','f4',('nch',),zlib=zlib,fill_value=MISSING)
            mr.long_name = 'refractive index of water' 
            mr.missing_value = np.float32(MISSING)
            mr.units         = "None"
            mr[:]            = self.__dict__['mr']  

            u10m = nc.createVariable('U10m','f4',('nobs',),zlib=zlib,fill_value=MISSING)
            u10m.long_name = '10-meter wind speed U component' 
            u10m.missing_value = np.float32(MISSING)
            u10m.units         = "None"
            u10m[:]            = self.__dict__['U10m']   

            v10m = nc.createVariable('V10m','f4',('nobs',),zlib=zlib,fill_value=MISSING)
            v10m.long_name = '10-meter wind speed V component' 
            v10m.missing_value = np.float32(MISSING)
            v10m.units         = "None"
            v10m[:]            = self.__dict__['V10m']  

        if ('NOBM' in self.albedoType) and (self.sleave_adjust is True):
            sl = nc.createVariable('adjusted_sleave','f4',('nch','vza','vaa',),zlib=zlib,fill_value=MISSING)
            sl.long_name = 'transmittance adjusted sleave'
            sl.missing_value = np.float32(MISSING)
            sl.units         = 'None'
            sl[:]            = self.adjusted_sleave

        if self.aerosol:
            tau = nc.createVariable('TAU','f4',('lev','nch',),zlib=zlib,fill_value=MISSING)
            tau.long_name = 'layer aerosol extinction optical depth' 
            tau.missing_value = np.float32(MISSING)
            tau.units         = "None"
            tau[:]            = np.squeeze(self.tau)

            ssa = nc.createVariable('SSA','f4',('lev','nch',),zlib=zlib,fill_value=MISSING)
            ssa.long_name = 'layer aerosol single scattering albedo' 
            ssa.missing_value = np.float32(MISSING)
            ssa.units         = "None"
            ssa[:]            = np.squeeze(self.ssa)

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
    format   = 'NETCDF4_CLASSIC'

    rootDir  = '/discover/nobackup/projects/gmao/osse2/pub/c1440_NR/OBS/POLAR_LIDAR/CALIPSO/'
    if not os.path.exists(rootDir):
        rootDir = '/nobackup/3/pcastell/POLAR_LIDAR/CALIPSO/'

    outRoot = './self_benchmark_2p8p2/graspConfig_12_Osku_DrySU_V1_hyperspectral/'


    #albedoType   = 'BPDF'
    #albedoType = 'MODIS_BRDF_BPDF'
    #albedoType = 'LAMBERTIAN_BPDF'
    #albedoType = 'GissCX'
    #albedoType = 'CX'
    #albedoType = 'LAMBERTIAN'
    #albedoType = 'MODIS_BRDF'
    albedoType  = 'OCIGissCX_NOBM_CLOUD'
    
    aerosol        = True # true if you want aerosols in simulation
    simple_aerosol = True # true if exponential aerosol with constnat optical properties
    dark           = True  #use if you want lambertian surface with albedo = 0
    nstreams       = 20
    plane_parallel = True
    sleave_adjust  = True  # if surface leaving surface, do you want to adjust for the transmittance

    step = 10
    channels  = np.arange(350,710,step)  
    chdmin    = get_chd(channels.min())
    chdmax    = get_chd(channels.max())
    

    if (albedoType == 'LAMBERTIAN') and (dark is True):
        sfcname = 'nosurface'
    elif albedoType == 'MODIS_BRDF':
        sfcname = 'BRDF'
    elif ('NOBM' in albedoType) and (sleave_adjust is True):
        sfcname = albedoType + '_sleave_adjust'
    else:
        sfcname = albedoType

    if aerosol is False:
        aname = 'rayleigh'
    else:
        if simple_aerosol:
            #aname = 'rayleigh+simple_aerosol'
            aname = 'rayleigh+simple_aerosol'
        else:
            aname = 'rayleigh+aerosol'

    outDir    = '{}/benchmark_{}_{}'.format(outRoot,aname,sfcname)
    outFile   = '{}/calipso-g5nr.vlidort.vector.{}.{}_{}_{}.nc4'.format(outDir,albedoType,chdmin,chdmax,step)
        
    rcFile   = 'rc/Aod_EOS.{}.rc'
    verbose  = True

    # Initialize VLIDORT class getting aerosol optical properties
    # -----------------------------------------------------------
    vlidort = BENCHMARK(albedoType,channels,outFile,rootDir,rcFile,
                             verbose=verbose,aerosol=aerosol,dark=dark,
                             simple_aerosol=simple_aerosol,
                             nstreams=nstreams,
                             plane_parallel=plane_parallel,
                             sleave_adjust=sleave_adjust)

       
    vlidort.runVLIDORT()
