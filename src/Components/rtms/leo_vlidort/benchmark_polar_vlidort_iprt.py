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


rootDir = '/home/pcastell/workspace/vlidort_benchmark2p8p2'

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

MISSING = -1.e+20

class BENCHMARK(POLAR_VLIDORT):
    """
    Everything needed for calling VLIDORT
    """
    def __init__(self,case,outDir,verbose=True):

        self.rootDir = rootDir
        self.outDir  = outDir
        self.verbose = verbose

        # get filesnames
        self.getFilenames(case)
            
        # set up cases
        if case == 'b3':
            self.setupB3()
        if case == 'brdf1':
            self.setupBRDF1()
        if case == 'brdf_bpdf1':
            self.setupBRDF_BPDF1()
        if case == 'gisscx':
            self.setupGISSCX()

        # calc angles
        self.calcAngles()

 
    # ---
    def getFilenames(self,case):
        self.inDir           = '{}/{}_dat/'.format(self.rootDir,case)
        self.outFile    = '{}/{}_benchmark.nc4'.format(self.outDir,case.upper())

    # ---
    def setupB3(self):
        self.SZA = 30.0
        self.SAA = 0.0 
        nvza = 18
        nvaa = 73
        self.VZA = np.arange(nvza)*5
        self.VAA = np.arange(nvaa)*5

        self.depol_ratio = 0.03
        self.nMom = 1000
        self.nPol = 6
        self.km = 30
        self.nstreams = 20
        self.plane_parallel = True
        self.albedoType = 'LAMBERTIAN'
        self.albedo = 0.0
        self.channel = 350.0

        # Read in ROT
        ze = []
        ROT = []
        f = open('{}/tau_rayleigh_350.dat'.format(self.inDir))
        f.readline() #header
        f.readline() #header
        for l in f:
            k,tau = l.split()
            ROT.append(float(tau))
            ze.append(float(k)+1)

        f.close()

        ze.append(0)

        self.ze = np.array(ze)
        self.ROT = np.array(ROT)
        
        self.ROT.shape = (self.km,1,1) #(km,nobs,nch)
        self.ze.shape = (self.km+1,1)

        # read molecular absorption
        MAOT = []
        f = open('{}/tau_molabs_350.dat'.format(self.inDir))
        f.readline() #header
        f.readline() #header
        for l in f:
            k, tau = l.split()
            MAOT.append(float(tau))

        f.close()

        self.alpha = np.array(MAOT)
        self.alpha.shape = (self.km,1,1) #(km,nobs,nch)

        # read AOT
        AOT = []
        f = open('{}/tau_aerosol.dat'.format(self.inDir))
        f.readline() #header
        f.readline() #header
        for l in f:
            k, tau = l.split()
            AOT.append(float(tau))

        f.close()

        self.tau = np.array(AOT)
        self.tau.shape = (self.km,1,1) #(km,nch,nobs)

        # read aerosol SSA
        nc = Dataset('{}/sizedistr_spheroid.cdf'.format(self.inDir))
        ssa = nc.variables['ssa'][:]
        nc.close()

        self.ssa = np.repeat(ssa,self.km)
        self.ssa.shape = (self.km,1,1) #(km,nch,nobs)

        # read expansion coefficients
        pmom = np.zeros([self.km,1,1,self.nMom,self.nPol])
        f = open('{}/xk_a4.txt'.format(self.inDir))
        f.readline() #header
        f.readline() #header
        for l in f:
            a, a1k, a2k, a3k, a4k, b1k, b2k = l.split()
            a = int(a)

            pmom[:,0,0,a,0] = float(a1k)
            pmom[:,0,0,a,1] = -1.*float(b1k)
            pmom[:,0,0,a,2] = float(a3k)
            pmom[:,0,0,a,3] = -1.*float(b2k)
            pmom[:,0,0,a,4] = float(a2k)
            pmom[:,0,0,a,5] = float(a4k)

        f.close()

        self.pmom = pmom  #(km,nch,nobs,nMom,nPol)

    # ---
    def setupBRDF1(self):
        self.SZA = 45.0
        self.SAA = 0.0
        nvza = 17
        nvaa = 73
        self.VZA = np.arange(nvza)*5
        self.VAA = np.arange(nvaa)*5

        self.depol_ratio = np.array([0.03])
        self.ROT = np.array([0.1])
        self.ze  = np.array([1,0])
        self.nMom = 32
        self.nPol = 6
        self.km = 1
        self.nstreams = 16
        self.plane_parallel = True
        self.albedoType = 'MODIS_BRDF'
        # [nkernel,nch,nobs]
        self.kernel_wt = np.zeros([3,1,1])
        self.kernel_wt[0,:,:]= 0.17424166  #fiso
        self.kernel_wt[1,:,:]= 0.044916667 #fgeo
        self.kernel_wt[2,:,:]= 0.041216664 # fvol

        # [nparam, nch, nobs]
        self.RTLSparam = np.zeros([2,1,1])
        self.RTLSparam[0,:,:]    = 2
        self.RTLSparam[1,:,:]    = 1        
        self.channel = 350.0

        chs = str(int(self.channel))
        self.__dict__['Riso'+chs] = self.kernel_wt[0,0,0]
        self.__dict__['Rgeo'+chs] = self.kernel_wt[1,0,0]
        self.__dict__['Rvol'+chs] = self.kernel_wt[2,0,0]


        self.ROT.shape = (self.km,1,1) #(km,nobs,nch)
        self.ze.shape = (self.km+1,1)

        # molecular absorption
        self.alpha = np.zeros([1])
        self.alpha.shape = (self.km,1,1) #(km,nobs,nch)

        # aerosol optical depth
        self.tau = np.zeros([1])
        self.tau.shape = (self.km,1,1) #(km,nch,nobs)

        # aerosol SSA
        self.ssa = np.zeros([1])
        self.ssa.shape = (self.km,1,1) #(km,nch,nobs)

        # aerosol expansion coefficients
        #(km,nch,nobs,nMom,nPol)
        self.pmom = np.zeros([self.km,1,1,self.nMom,self.nPol])

    # ---
    def setupBRDF_BPDF1(self):
        self.SZA = 45.0
        self.SAA = 0.0
        nvza = 17
        nvaa = 73
        self.VZA = np.arange(nvza)*5
        self.VAA = np.arange(nvaa)*5

        self.depol_ratio = np.array([0.03])
        self.ROT = np.array([0.1])
        self.ze  = np.array([1,0])
        self.nMom = 32
        self.nPol = 6
        self.km = 1
        self.nstreams = 16
        self.plane_parallel = True
        self.albedoType = 'MODIS_BRDF_BPDF'
        # [nkernel,nch,nobs]
        self.kernel_wt = np.zeros([3,1,1])
        self.kernel_wt[0,:,:]= 0.17424166  #fiso
        self.kernel_wt[1,:,:]= 0.044916667 #fgeo
        self.kernel_wt[2,:,:]= 0.041216664 # fvol

        # [nparam, nch, nobs]
        self.RTLSparam = np.zeros([2,1,1])
        self.RTLSparam[0,:,:]    = 2
        self.RTLSparam[1,:,:]    = 1

        #BPDFparam(nparam,nch,nobs)
        self.BPDFparam = np.zeros([3,1,1])
        self.BPDFparam[0,0,:] = 1.5
        self.BPDFparam[1,0,:] = 0.78 #NDVI
        self.BPDFparam[2,0,:] = 6.86 #BPDFcoef

        self.channel = 350.0

        chs = str(int(self.channel))
        self.__dict__['Riso'+chs] = self.kernel_wt[0,0,0]
        self.__dict__['Rgeo'+chs] = self.kernel_wt[1,0,0]
        self.__dict__['Rvol'+chs] = self.kernel_wt[2,0,0]


        self.ROT.shape = (self.km,1,1) #(km,nobs,nch)
        self.ze.shape = (self.km+1,1)

        # molecular absorption
        self.alpha = np.zeros([1])
        self.alpha.shape = (self.km,1,1) #(km,nobs,nch)

        # aerosol optical depth
        self.tau = np.zeros([1])
        self.tau.shape = (self.km,1,1) #(km,nch,nobs)

        # aerosol SSA
        self.ssa = np.zeros([1])
        self.ssa.shape = (self.km,1,1) #(km,nch,nobs)

        # aerosol expansion coefficients
        #(km,nch,nobs,nMom,nPol)
        self.pmom = np.zeros([self.km,1,1,self.nMom,self.nPol])

    # ---
    def setupGISSCX(self):
        self.SZA = 45.0
        self.SAA = 0.0
        nvza = 17
        nvaa = 73
        self.VZA = np.arange(nvza)*5
        self.VAA = np.arange(nvaa)*5

        self.depol_ratio = np.array([0.03])
        self.ROT = np.array([0.1])
        self.ze  = np.array([1,0])
        self.nMom = 32
        self.nPol = 6
        self.km = 1
        self.nstreams = 16
        self.plane_parallel = True
        self.albedoType = 'OCIGissCX'
        # [nkernel,nch,nobs]
        self.kernel_wt = np.zeros([3,1,1])
        self.kernel_wt[0,:,:]= 0.17424166  #fiso
        self.kernel_wt[1,:,:]= 0.044916667 #fgeo
        self.kernel_wt[2,:,:]= 0.041216664 # fvol

        self.mr = np.array([1.333])
        self.U10m = np.array([3.0])
        self.V10m = np.array([4.0])

        self.channel = 350.0

        self.ROT.shape = (self.km,1,1) #(km,nobs,nch)
        self.ze.shape = (self.km+1,1)

        # molecular absorption
        self.alpha = np.zeros([1])
        self.alpha.shape = (self.km,1,1) #(km,nobs,nch)

        # aerosol optical depth
        self.tau = np.zeros([1])
        self.tau.shape = (self.km,1,1) #(km,nch,nobs)

        # aerosol SSA
        self.ssa = np.zeros([1])
        self.ssa.shape = (self.km,1,1) #(km,nch,nobs)

        # aerosol expansion coefficients
        #(km,nch,nobs,nMom,nPol)
        self.pmom = np.zeros([self.km,1,1,self.nMom,self.nPol])


    # ---
    def computeAtmos(self):

        pe, ze, te = getEdgeVars(self)

        self.pe = pe # (km,nobs)
        self.ze = ze 
        self.te = te

    # --
    def calcAngles(self):
        raa = self.VAA - self.SAA
        I = raa < 0.
        raa[I] = 360.0 + raa[I]

        self.RAA = raa

    def runVLIDORT(self):
        """
        Calls VLIDORT 
        """
        flux_factor = np.ones([1,1]) #solar flux F0 (nch,nobs)
        pe   = self.ze*0.0
        ze   = self.ze*1e3
        te   = self.ze*0.0

        nlev = self.km
        tau  = self.tau   
        ssa  = self.ssa   
        pmom = self.pmom 
        alpha = self.alpha
        

        nvza = len(self.VZA)
        nraa = len(self.RAA)
        sza  = self.SZA

        self.I = np.ones([nvza,nraa])*MISSING
        self.Q = np.ones([nvza,nraa])*MISSING
        self.U = np.ones([nvza,nraa])*MISSING
        self.reflectance = np.ones([nvza,nraa])*MISSING
        self.surf_reflectance = np.ones([nvza,nraa])*MISSING
        self.BR_Q = np.ones([nvza,nraa])*MISSING
        self.BR_U = np.ones([nvza,nraa])*MISSING

        ROT = self.ROT
        depol_ratio = self.depol_ratio 
        
        # Get VLIDORT wrapper name from dictionary
        vlidortWrapper = WrapperFuncs[self.albedoType]

        # loop through viewing angles and run VLIDORT
        for ivza,vza in enumerate(self.VZA[0:1]):
            for iraa,raa in enumerate(self.RAA[0:1]):
                    
                # get args list for each surface model
                if self.albedoType == 'MODIS_BRDF':
                    kernel_wt = self.kernel_wt[:,:,0]
                    param     = self.RTLSparam[:,:,0]                
                    

                    args = [self.channel, self.nstreams, self.plane_parallel, ROT, depol_ratio, tau, ssa, pmom, 
                            pe, ze, te, 
                            kernel_wt, param, 
                            sza, raa, vza, 
                            MISSING,
                            self.verbose]

                    # Call VLIDORT wrapper function
                    I, reflectance, surf_reflectance, Q, U, BR_Q, BR_U, rc = vlidortWrapper(*args)      
                    print('I',I)
                    print('rc',rc)                    

                elif self.albedoType == 'MODIS_BRDF_BPDF':
                    kernel_wt = self.kernel_wt[:,:,0]
                    param     = self.RTLSparam[:,:,0]
                    bpdf      = self.BPDFparam[:,:,0]

                    args = [self.channel, self.nstreams, self.plane_parallel, ROT, depol_ratio, tau, ssa, pmom,
                            pe, ze, te,
                            kernel_wt, param, bpdf,
                            sza, raa, vza,
                            MISSING,
                            self.verbose]

                    # Call VLIDORT wrapper function
                    I, reflectance, surf_reflectance, Q, U, BR_Q, BR_U, rc = vlidortWrapper(*args)
                    print('I',I)
                    print('rc',rc)

                elif ('CX' in self.albedoType) and (self.albedoType != 'OCIGissCX_NOBM_CLOUD'):
                    args = [self.channel, self.nstreams, self.plane_parallel, ROT, depol_ratio, alpha, tau, ssa, pmom,
                            pe, ze, te,
                            self.U10m, self.V10m, self.mr,
                            sza, raa, vza,
                            flux_factor,
                            MISSING,
                            self.verbose]                    

                    # Call VLIDORT wrapper function
                    I, reflectance, surf_reflectance, Q, U, BR_Q, BR_U, rc = vlidortWrapper(*args)
                    print('I',I)
                    print('rc',rc)                    
                elif self.albedoType == 'LAMBERTIAN':
                    albedo = self.albedo
                    
                    args = [self.channel, self.nstreams, self.plane_parallel, ROT, depol_ratio, alpha, tau, ssa, pmom, 
                            pe, ze, te, 
                            albedo, 
                            sza, raa, vza,
                            flux_factor,
                            MISSING,
                            self.verbose]

                    # Call VLIDORT wrapper function
                    I, reflectance, Q, U, rc = vlidortWrapper(*args)  
                    surf_reflectance = albedo
                    print('I',I)
                    print('rc',rc)

                self.I[ivza,iraa] = np.squeeze(I)
                self.reflectance[ivza,iraa] = np.squeeze(reflectance)
                self.surf_reflectance[ivza,iraa] = np.squeeze(surf_reflectance)
                self.Q[ivza,iraa] = np.squeeze(Q)
                self.U[ivza,iraa] = np.squeeze(U) 
                if self.albedoType != 'LAMBERTIAN':
                    self.BR_Q[ivza,iraa] = np.squeeze(BR_Q)
                    self.BR_U[ivza,iraa] = np.squeeze(BR_U) 

        
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
        nc.SZA    = self.SZA
        nc.SAA    = self.SAA

        if 'BRDF' in self.albedoType:
            nc.Riso = self.kernel_wt[0,0,0]
            nc.Rgeo = self.kernel_wt[1,0,0]
            nc.Rvol = self.kernel_wt[2,0,0]
     
        # Create dimensions
        # -----------------
        nt = nc.createDimension('vza',len(self.VZA))
        ls = nc.createDimension('vaa',len(self.VAA))
        nz = nc.createDimension('lev',self.km)
        nze = nc.createDimension('leve',self.km+1)
        no = nc.createDimension('nobs',1)   
        # if self.aerosol:
        #     nr = nc.createDimension('radius',len(self.R))  
        #     na = nc.createDimension('angle',len(self.angle))    


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
        ref = nc.createVariable('toa_reflectance','f4',('vza','vaa',),zlib=zlib,fill_value=MISSING)
        ref.standard_name = '%.2f nm TOA Reflectance' %self.channel
        ref.long_name     = '%.2f nm reflectance at the top of the atmosphere' %self.channel
        ref.missing_value = np.float32(MISSING)
        ref.units         = "None"
        ref[:]            = self.reflectance

        i = nc.createVariable('I','f4',('vza','vaa',),zlib=zlib,fill_value=MISSING)
        i.standard_name = '%.2f nm TOA I' %self.channel
        i.long_name     = '%.2f nm intensity at the top of the atmosphere' %self.channel
        i.missing_value = np.float32(MISSING)
        i.units         = "W m-2 sr-1 nm-1"
        i[:]            = self.I

        q = nc.createVariable('Q','f4',('vza','vaa',),zlib=zlib,fill_value=MISSING)
        q.standard_name = '%.2f nm TOA Q' %self.channel
        q.long_name     = '%.2f nm Q-component of the stokes vector at the top of the atmopshere' %self.channel
        q.missing_value = np.float32(MISSING)
        q.units         = "W m-2 sr-1 nm-1"
        q[:]            = self.Q    

        u = nc.createVariable('U','f4',('vza','vaa',),zlib=zlib,fill_value=MISSING)
        u.standard_name = '%.2f nm TOA U' %self.channel
        u.long_name     = '%.2f nm U-component of the stokes vector at the top of the atmopshere' %self.channel
        u.missing_value = np.float32(MISSING)
        u.units         = "W m-2 sr-1 nm-1"
        u[:]            = self.U

        sref = nc.createVariable('surf_reflectance','f4',('vza','vaa',),zlib=zlib,fill_value=MISSING)
        sref.standard_name = '%.2f nm Surface Reflectance' %self.channel
        sref.long_name     = '%.2f nm Bi-Directional Surface Reflectance' %self.channel
        sref.missing_value = np.float32(MISSING)
        sref.units         = "None"
        sref[:]            = self.surf_reflectance

        sref = nc.createVariable('surf_reflectance_Q','f4',('vza','vaa',),zlib=zlib,fill_value=MISSING)
        sref.standard_name = '%.2f nm Surface Reflectance Q' %self.channel
        sref.long_name     = '%.2f nm Bi-Directional Surface Reflectance Q' %self.channel
        sref.missing_value = np.float32(MISSING)
        sref.units         = "None"
        sref[:]            = self.BR_Q

        sref = nc.createVariable('surf_reflectance_U','f4',('vza','vaa',),zlib=zlib,fill_value=MISSING)
        sref.standard_name = '%.2f nm Surface Reflectance U' %self.channel
        sref.long_name     = '%.2f nm Bi-Directional Surface Reflectance U' %self.channel
        sref.missing_value = np.float32(MISSING)
        sref.units         = "None"
        sref[:]            = self.BR_U


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

    outRoot = './self_benchmark_2p8p2/'
    case    = 'b3'
    outDir    = '{}/{}'.format(outRoot,case)
        
    # Initialize VLIDORT class getting aerosol optical properties
    # -----------------------------------------------------------
    vlidort = BENCHMARK(case,outDir)

       
    vlidort.runVLIDORT()
    vlidort.writeNC()
