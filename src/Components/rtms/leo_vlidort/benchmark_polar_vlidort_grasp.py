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

GRASP_phase = 'Expansion_Coefs/GRASPe438668_correctNorm_simpleAerosol_PMonly.txt'

nMom     = 300
#nMom     = 30 # GRASP Phase function
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
MieVarsNames = ['ext','scatext','backscat','aback_sfc','aback_toa','depol','ext2back','tau','ssa','g']
MieVarsUnits = ['km-1','km-1','km-1 sr-1','sr-1','sr-1','unitless','sr','unitless','unitless','unitless']

META    = ['DELP','PS','RH','AIRDENS','trjLon','trjLat','isotime']
AERNAMES = VNAMES_SU #+ VNAMES_SS + VNAMES_OC + VNAMES_BC + VNAMES_DU
#AERNAMES = ['DU001']
SDS_AER = META + AERNAMES
MISSING = -1.e+20

class BENCHMARK(POLAR_VLIDORT):
    """
    Everything needed for calling VLIDORT
    """
    def __init__(self,albedoType,channel,outFile,rootDir,rcFile,
                verbose=False,aerosol=True,dark=False,simple_aerosol=False,
                emde=False,grasp=False,grasp_higher=False,rayleigh=False,
                nstreams=12,plane_parallel=True):
        self.outFile = outFile
        self.albedoType = albedoType
        self.channel = channel
        self.verbose = verbose
        self.nMom    = nMom
        self.nPol    = nPol
        self.nobs    = 1
        self.nch     = 1
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

        # get filesnames
        self.getFilenames()
            
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

        # Read Model Data
        self.readSampledGEOS()
        # Calculate atmospheric profile properties needed for Rayleigh calc      
        self.computeAtmos()
        # Calculate aerosol optical properties
        # set to zero if aerosol is false
        self.computeMie()          
        # if self.aerosol:
        #     # calculate size distribution
        #     self.sizeDistribution()

        # Calculate Scene Geometry
        self.VZA = VZA
        self.VAA = VAA
        self.SZA = np.array(SZA)
        self.SAA = np.array(SAA)
        self.calcAngles()

    def getFilenames(self):
        inDir           = '{}/LevelB/Y2006/M08/'.format(self.rootDir)
        self.inFile     = '{}/calipso-g5nr.lb2.aer_Nv.20060801_00z.nc4'.format(inDir)
        inDir           = '{}/LevelB/surface/BRDF/MCD43C1/006/Y2006/M08'.format(self.rootDir)
        self.brdfFile   = '{}/calipso-g5nr.lb2.brdf.20060801_00z.nc4'.format(inDir)
        inDir           = '{}/LevelB/surface/SurfLER/Y2006/M08'.format(self.rootDir)
        self.lerFile    = '{}/calipso-g5nr.lb2.omi_ler.20060801_00z.nc4'.format(inDir)

    #---
    def computeMie(self):
        """
        Computes aerosol optical quantities 
        """
        if self.simple_aerosol and self.aerosol:
            self.RH[:] = 0.00

            tauO,ssa,g,pmom = getAOPvector(self,self.channel,
                                 vnames=self.AERNAMES,
                                 Verbose=True,
                                 rcfile=self.rcFile,
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

            if self.emde:
                # read in Emde et al. legendre coefficients
                f = open('/home/pcastell/workspace/vlidort_benchmark/a3_dat/xk_a3.txt','r')
                f.readline() #header
                f.readline() #header

                pmom[:] = 0.0
                ich = 0
                iob = 0
                n = 0
                for l in f:
                    k,a1k,a2k,a3k,a4k,b1k,b2k = l.split()                    
                    pmom[:,ich,iob,n,0] = float(a1k)
                    pmom[:,ich,iob,n,1] = -1.0*float(b1k)
                    pmom[:,ich,iob,n,2] = float(a3k)
                    pmom[:,ich,iob,n,3] = -1.0*float(b2k)
                    pmom[:,ich,iob,n,4] = float(a2k)
                    pmom[:,ich,iob,n,5] = float(a4k)
                    n = n + 1
            elif self.grasp:
                # Read in phase function generated by GRASP
                # fit legendre coefficients
                f = open('PM_GRASP_simple_aerosol.txt','r')
                # headers
                f.readline()
                f.readline()
                f.readline()
                f.readline()
                f.readline()
                
                pmom[:] = 0.0
                ich = 0
                iob = 0
                n = 0       
                P11, P12, P33, P34, ANG = np.zeros(181), np.zeros(181),np.zeros(181),np.zeros(181),np.zeros(181)   
                for l in f:
                    k,ang,p11,p12,p22,p33,p34,p44 = l.split()
                    P11[n] = float(p11)
                    P12[n] = float(p12)
                    P33[n] = float(p33)
                    P34[n] = float(p34)
                    ANG[n] = float(ang)

                    n = n+ 1

                mu    = np.cos(np.radians(ANG))
                deg   = 20
                coef11  = np.polynomial.legendre.legfit(mu,P11,deg)
                coef12  = np.polynomial.legendre.legfit(mu,P12,deg)
                coef33  = np.polynomial.legendre.legfit(mu,P33,deg)
                coef34  = np.polynomial.legendre.legfit(mu,P34,deg)
                coef11[0] = 1.0

                for n in range(deg+1):
                    pmom[:,ich,iob,n,0] = coef11[n]

                    pmom[:,ich,iob,n,2] = coef33[n]
                    pmom[:,ich,iob,n,3] = coef34[n]
                    pmom[:,ich,iob,n,4] = coef11[n]
                    pmom[:,ich,iob,n,5] = coef33[n]

                for n in range(deg+1):                    
                    pmom[:,ich,iob,n,1] = coef12[n]

                self.coef11 = coef11
                self.coef12 = coef12
                self.coef33 = coef33
                self.coef34 = coef34

            elif self.grasp_higher:
                # Read in coefficients calculated with Mischenko's code
                f = open('{}.expan_coeff'.format(GRASP_phase),'r')
                f.readline() #header
                pmom[:] = 0
                ich = 0
                iob = 0
                for l in f:
                    n, P11, P22, P33, P44, P12, P34 = l.split()
                    n = int(n)
                    pmom[:,ich,iob,n,0] = float(P11)
                    pmom[:,ich,iob,n,1] = -1.*float(P12)
                    pmom[:,ich,iob,n,2] = float(P33)
                    pmom[:,ich,iob,n,3] = float(P34)
                    pmom[:,ich,iob,n,4] = float(P22)
                    pmom[:,ich,iob,n,5] = float(P44)

                f.close()


                # Read in reconstructed function calculated with Mischenko's code
                f = open('{}.expan_matr'.format(GRASP_phase),'r')
                p11 = []
                p22 = []
                p33 = []
                p44 = []
                p12 = []
                p34 = []
                angle = []
                for l in f:
                    n, P11, P22, P33, P44, P12, P34 = l.split()
                    angle.append(float(n))
                    p11.append(float(P11))
                    p22.append(float(P22))
                    p33.append(float(P33))
                    p44.append(float(P44))
                    p12.append(float(P12))
                    p34.append(float(P34))

                self.angle = np.array(angle)
                self.p11 = np.array(p11)
                self.p22 = np.array(p22)
                self.p33 = np.array(p33)
                self.p44 = np.array(p44)
                self.p12 = np.array(p12)
                self.p34 = np.array(p34)

                # Read in phase function generated by GRASP
                f = open(GRASP_phase,'r')
                
                n = 0       
                P11, P12, P22, P33, P34, P44, ANG = np.zeros(181), np.zeros(181),np.zeros(181),np.zeros(181),np.zeros(181),np.zeros(181),np.zeros(181)     
                for l in f:
                    ang,p11,p12,p22,p33,p34,p44 = l.split()
                    P11[n] = float(p11)
                    P22[n] = float(p22)
                    P33[n] = float(p33)                    
                    P44[n] = float(p44)
                    P12[n] = float(p12)
                    P34[n] = float(p34)
                    ANG[n] = float(ang)

                    n = n+ 1

                self.P11 = P11
                self.P22 = P22
                self.P33 = P33
                self.P44 = P44
                self.P12 = P12
                self.P34 = P34

            elif self.rayleigh:
                pmom[:] = 0
                ich = 0
                iob = 0
                coef11 = [1.0, 0.0, 0.4778325]
                coef12 = [0.0, 0.0, -1.170446]
                coef22 = [0.0, 0.0, 2.866995]
                coef33 = [0.0, 0.0, 0.0]
                coef34 = [0.0, 0.0, 0.0]
                coef44 = [0.0, 1.38916256, 0.0]
                for n in range(3):
                    pmom[:,ich,iob,n,0] = coef11[n]
                    pmom[:,ich,iob,n,1] = coef12[n]
                    pmom[:,ich,iob,n,2] = coef33[n]
                    pmom[:,ich,iob,n,3] = coef34[n]
                    pmom[:,ich,iob,n,4] = coef22[n]
                    pmom[:,ich,iob,n,5] = coef44[n]


        elif self.aerosol:
            tau,ssa,g,pmom = getAOPvector(self,self.channel,
                                 vnames=self.AERNAMES,
                                 Verbose=True,
                                 rcfile=self.rcFile,
                                 nMom=self.nMom)
        else:
            tau = np.zeros([self.km,self.nch,self.nobs])
            ssa = np.zeros([self.km,self.nch,self.nobs])
            g   = np.zeros([self.km,self.nch,self.nobs])
            pmom = np.zeros([self.km,self.nch,self.nobs,self.nMom,self.nPol])

        self.tau = tau  #(km,nch,nobs)
        self.ssa = ssa  #(km,nch,nobs)
        self.g   = g    #(km,nch,nobs)
        self.pmom = pmom  #(km,nch,nobs,nMom,nPol)

        # if self.aerosol:
            # vol, area, refr, refi, reff = getAOPint(self,self.channel,
            #                                         vnames=self.AERNAMES,
            #                                         Verbose=True,
            #                                         rcfile=self.rcFile)
            # self.vol  = vol
            # self.area = area
            # self.refr = refr
            # self.refi = refi
            # self.reff = reff

            # if not self.grasp_higher:
            #     tau,ssa,g,pmom = getAOPvector(self,self.channel,
            #                      vnames=self.AERNAMES,
            #                      Verbose=True,
            #                      rcfile='rc/Aod_EOS.800.old.rc',
            #                      nMom=300)                
            #     # calculate phase function from legendre coefficients
            #     angle = np.arange(181)
            #     mu    = np.cos(np.radians(angle))
            #     m = mu.shape

            #     p11, p12, p33, p34 = (np.zeros(m),np.zeros(m),np.zeros(m),np.zeros(m))
            #     k = 0
            #     ich = 0
            #     iob = 0
            #     #for n in range(self.nMom):
            #     for n in range(300):
            #         print "Adding moment ", n
            #         P = legendre(n)
                        
            #         p11 += pmom[k,ich,iob,n,0]* P(mu)
            #         p12 += pmom[k,ich,iob,n,1] * P(mu)
            #         p33 += pmom[k,ich,iob,n,2] * P(mu)
            #         p34 += pmom[k,ich,iob,n,3] * P(mu)

            #     self.angle = angle
            #     self.p11 = p11
            #     self.p12 = p12
            #     self.p33 = p33
            #     self.p34 = p34
                  
    # --
    def readSampledWindCX(self):
        """
        Cox-munk inputs
        """
        self.mr = np.array([1.333])
        self.U10m = np.array([3.0])
        self.V10m = np.array([4.0])


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

    def readSampledLER(self):

        if self.dark:
            self.albedo = np.array([0.0])  
        else:
            self.albedo = np.array([0.3])  

        self.albedo.shape = (1,1)


    def sizeDistribution(self):
        """ 
        Get size aerosol size distribution from 
        GEOS-5 aerosol mixing ratios

        Based on P. Colarco's optics table calculations.
        Note: hard coded for these versions of the optics tables.
        May not work for other versions.

        filename_optical_properties_DU: ExtData/optics_DU.v15_5.nc
        filename_optical_properties_SS: ExtData/optics_SS.v3_5.nc
        filename_optical_properties_BC: ExtData/optics_BC.v1_5.nc
        filename_optical_properties_OC: ExtData/optics_OC.v1_5.nc
        filename_optical_properties_SU: ExtData/optics_SU.v1_5.nc
        """

        # Create master bins for all levels
        # ---------------------------------
        rmin      = 0.001e-6   #meters 
        rmax      = self.getRMAX()
        R, DR, RLOW, RUP = self.logBins(0.5*rmin,1.1*rmax)
        self.R = R
        self.DR = DR
        self.RLOW = RLOW
        self.RUP  = RUP


        # Lognormal Species: BC, OC, and SU
        # ---------------------------------
        
        # BC
        r0        = 0.0118e-6  #meters
        rmin      = 0.001e-6   #meters 
        rmax0     = 5.0e-6     #meters
        sigma     = 2.00
        rhop0     = 1000       # Density of dry particles [kg m-3]

        spc = 'BCPHOBIC'
        if spc in self.AERNAMES:
            self.logNormalDistribution(spc,r0,rmin,rmax0,sigma,rhop0)
        spc = 'BCPHILIC'
        if spc in self.AERNAMES:
            self.logNormalDistribution(spc,r0,rmin,rmax0,sigma,rhop0)

        # OC
        r0        = 0.0212e-6  #meters
        rmin      = 0.005e-6   #meters 
        rmax0     = 5.0e-6     #meters
        sigma     = 2.12
        rhop0     = 1800       #[kg m-3]

        spc = 'OCPHOBIC'
        if spc in self.AERNAMES:
            self.logNormalDistribution(spc,r0,rmin,rmax0,sigma,rhop0)
        spc = 'OCPHILIC'
        if spc in self.AERNAMES:
            self.logNormalDistribution(spc,r0,rmin,rmax0,sigma,rhop0)

        # SU
        r0        = 0.0695e-6  #meters
        rmin      = 0.005e-6   #meters 
        rmax0     = 5.0e-6     #meters
        sigma     = 1.77
        rhop0     = 1700       #[kg m-3]

        spc = 'SU'
        if 'SO4' in self.AERNAMES:
            self.logNormalDistribution(spc,r0,rmin,rmax0,sigma,rhop0)

        # Dust
        # ----------
        if 'DU001' in self.AERNAMES:
            self.dustDistribution()

        # Sea Salt
        # ------------
        if 'SS001' in self.AERNAMES:
            self.seasaltDistribution()

        # Mixture Distribution
        # -------------
        #self.TOTdist = self.BCPHILICdist + self.BCPHOBICdist + self.OCPHILICdist + self.OCPHOBICdist + self.SUdist + self.DUdist + self.SSdist
        #self.BCdist  = self.BCPHILICdist + self.BCPHOBICdist
        #self.OCdist  = self.OCPHILICdist + self.OCPHOBICdist
        distnames = ['BCPHILICdist','BCPHOBICdist','OCPHILICdist','OCPHOBICdist','SUdist','DUdist','SSdist']
        for spc in distnames:
            if hasattr(self,spc):
                if hasattr(self,'TOTdist'):
                    self.TOTdist += self.__dict__[spc]
                else:
                    self.TOTdist = self.__dict__[spc]

                if spc in ('BCPHILICdist','BCPHOBICdist'):
                    if hasattr(self,'BCdist'):
                        self.BCdist += self.__dict__[spc]
                    else:
                        self.BCdist = self.__dict__[spc]

                if spc in ('OCPHILICdist','OCPHOBICdist'):
                    if hasattr(self,'BCdist'):
                        self.OCdist += self.__dict__[spc]
                    else:
                        self.OCdist = self.__dict__[spc]        

        # Column size distribution
        self.colTOTdist = np.zeros([self.nobs,len(self.R)])
        dz = self.ze[:-1,0] - self.ze[1:,0]  
        for t in range(self.nobs):
            for r in range(len(self.R)):
                self.colTOTdist[t,r] = np.sum(self.TOTdist[t,:,r]*dz)
        


        # convert to dV/dlnR [microns^3/microns^2] 
        distnames = ['BCdist','OCdist','SUdist','DUdist','SSdist','TOTdist']
        nlev      = 72
        for spc in distnames:
            if hasattr(self,spc):  
                temp = np.zeros(self.__dict__[spc].shape)  
                for t in range(self.nobs):
                    for k in range(nlev):
                        temp[t,k,:] = self.__dict__[spc][t,k,:]*self.R*1e6

                self.__dict__[spc] = temp

        spc = 'colTOTdist'
        temp = np.zeros(self.__dict__[spc].shape)  
        for t in range(self.nobs):
            temp[t,:] = self.__dict__[spc][t,:]*self.R*1e6

        self.__dict__[spc] = temp

    def rot2scatplane(self):
        """
        Rotate meridonal plane (VLIDORT) outputs to scattering plane (GRASP)
        adapted from Kirk's code rot2scatplane.pro
        """
        
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
        self.Q_out = np.zeros([nvza,nraa])
        self.U_out = np.zeros([nvza,nraa])
        self.Qsurf_out = np.zeros([nvza,nraa])
        self.Usurf_out = np.zeros([nvza,nraa])        
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
                q_in = self.Q[ivza,iraa]
                u_in = self.U[ivza,iraa]
                q_out= q_in*cos2i2 - u_in*sin2i2
                u_out= q_in*sin2i2 + u_in*cos2i2

                self.Q_out[ivza,iraa] = -1.*q_out
                self.U_out[ivza,iraa] = -1.*u_out

                q_in = self.BR_Q[ivza,iraa]
                u_in = self.BR_U[ivza,iraa]
                q_out= q_in*cos2i2 - u_in*sin2i2
                u_out= q_in*sin2i2 + u_in*cos2i2                
                self.Qsurf_out[ivza,iraa] = -1.*q_out
                self.Usurf_out[ivza,iraa] = -1.*u_out



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

        # Calculate ROT
        args = [self.channel, pe, ze, te, MISSING, self.verbose]
        vlidortWrapper = WrapperFuncs['ROT_CALC']
        ROT, depol_ratio, rc = vlidortWrapper(*args) 
        depol_ratio = np.array([0.03])
        

        # hack to scale to GRASP
        rod = 1.72484312e-02        
        ROT = ROT*rod/np.sum(ROT)
        self.ROT = np.squeeze(ROT).T

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
                    

                    args = [self.channel, self.nstreams, self.plane_parallel, ROT, depol_ratio, tau, ssa, pmom, 
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

                    args = [self.channel,
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

                    args = [self.channel,
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

                    args = [self.channel,
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
                
                elif 'CX' in self.albedoType:
                    args = [self.channel, self.nstreams, self.plane_parallel, ROT, depol_ratio, tau, ssa, pmom,
                            pe, ze, te,
                            self.U10m, self.V10m, self.mr,
                            sza, raa, vza,
                            MISSING,
                            self.verbose]

                    I, reflectance, surf_reflectance, Q, U, BR_Q, BR_U, rc = vlidortWrapper(*args)
                elif self.albedoType == 'LAMBERTIAN':
                    albedo = self.albedo
                    
                    args = [self.channel, self.nstreams, self.plane_parallel, ROT, depol_ratio, tau, ssa, pmom, 
                            pe, ze, te, 
                            albedo, 
                            sza, raa, vza, 
                            MISSING,
                            self.verbose]

                    # Call VLIDORT wrapper function
                    I, reflectance, Q, U, rc = vlidortWrapper(*args)  
                    surf_reflectance = albedo

                self.I[ivza,iraa] = np.squeeze(I)
                self.reflectance[ivza,iraa] = np.squeeze(reflectance)
                self.surf_reflectance[ivza,iraa] = np.squeeze(surf_reflectance)
                self.Q[ivza,iraa] = np.squeeze(Q)
                self.U[ivza,iraa] = np.squeeze(U) 
                if self.albedoType != 'LAMBERTIAN':
                    self.BR_Q[ivza,iraa] = np.squeeze(BR_Q)
                    self.BR_U[ivza,iraa] = np.squeeze(BR_U) 

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
        nt = nc.createDimension('vza',len(self.VZA))
        ls = nc.createDimension('vaa',len(self.VAA))
        nz = nc.createDimension('lev',km)
        nze = nc.createDimension('leve',km+1)
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
        vaa.missing_value = MISSING
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

        q = nc.createVariable('Q_scatplane','f4',('vza','vaa',),zlib=zlib,fill_value=MISSING)
        q.standard_name = '%.2f nm TOA Q rotated to scattering plane' %self.channel
        q.long_name     = '%.2f nm Q-component of the stokes vector at the top of the atmopshere' %self.channel
        q.missing_value = np.float32(MISSING)
        q.units         = "W m-2 sr-1 nm-1"
        q[:]            = self.Q_out

        u = nc.createVariable('U_scatplane','f4',('vza','vaa',),zlib=zlib,fill_value=MISSING)
        u.standard_name = '%.2f nm TOA U rotated to scattering plane' %self.channel
        u.long_name     = '%.2f nm U-component of the stokes vector at the top of the atmopshere' %self.channel
        u.missing_value = np.float32(MISSING)
        u.units         = "W m-2 sr-1 nm-1"
        u[:]            = self.U_out
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

        sref = nc.createVariable('surf_reflectance_Q_scatplane','f4',('vza','vaa',),zlib=zlib,fill_value=MISSING)
        sref.standard_name = '%.2f nm Surface Reflectance Q' %self.channel
        sref.long_name     = '%.2f nm Bi-Directional Surface Reflectance Q' %self.channel
        sref.missing_value = np.float32(MISSING)
        sref.units         = "None"
        sref[:]            = self.Qsurf_out

        sref = nc.createVariable('surf_reflectance_U_scatplane','f4',('vza','vaa',),zlib=zlib,fill_value=MISSING)
        sref.standard_name = '%.2f nm Surface Reflectance U' %self.channel
        sref.long_name     = '%.2f nm Bi-Directional Surface Reflectance U' %self.channel
        sref.missing_value = np.float32(MISSING)
        sref.units         = "None"
        sref[:]            = self.Usurf_out


        rot = nc.createVariable('ROT','f4',('lev',),zlib=zlib,fill_value=MISSING)
        rot.long_name = '%.2f nm Rayleigh Optical Thickness' %self.channel
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
            mr = nc.createVariable('mr','f4',('nobs',),zlib=zlib,fill_value=MISSING)
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

        if self.aerosol:
            tau = nc.createVariable('TAU','f4',('lev',),zlib=zlib,fill_value=MISSING)
            tau.long_name = 'layer aerosol extinction optical depth' 
            tau.missing_value = np.float32(MISSING)
            tau.units         = "None"
            tau[:]            = np.squeeze(self.tau)

            ssa = nc.createVariable('SSA','f4',('lev',),zlib=zlib,fill_value=MISSING)
            ssa.long_name = 'layer aerosol single scattering albedo' 
            ssa.missing_value = np.float32(MISSING)
            ssa.units         = "None"
            ssa[:]            = np.squeeze(self.ssa)

            # refi = nc.createVariable('REFI','f4',('lev',),zlib=zlib,fill_value=MISSING)
            # refi.long_name = 'layer aerosol imaginary refractive index' 
            # refi.missing_value = np.float32(MISSING)
            # refi.units         = "None"
            # refi[:]            = np.squeeze(self.refi)

            # refr = nc.createVariable('REFR','f4',('lev',),zlib=zlib,fill_value=MISSING)
            # refr.long_name = 'layer aerosol real refractive index' 
            # refr.missing_value = np.float32(MISSING)
            # refr.units         = "None"
            # refr[:]            = np.squeeze(self.refr)

            # reff = nc.createVariable('REFF','f4',('lev',),zlib=False)
            # reff.long_name   = 'aerosol effective radius'
            # reff.units       = 'microns'
            # reff[:]          = np.squeeze(self.reff*1e6)

            # area = nc.createVariable('AREA','f4',('lev',),zlib=False)
            # area.long_name   = 'aerosol cross sectional area'
            # area.units       = 'm2 m-3'
            # area[:]          = np.squeeze(self.area)

            # vol = nc.createVariable('VOL','f4',('lev',),zlib=False)
            # vol.long_name   = 'aerosol volume'
            # vol.units       = 'm3 m-3'
            # vol[:]          = np.squeeze(self.vol)

            # rad = nc.createVariable('radius','f4',('radius',),zlib=False)
            # rad.short_name        = 'R'
            # rad.long_name   = 'aerosol size distribution radius'
            # rad.units       = 'microns'
            # rad[:]          = self.R*1e6

            # drad = nc.createVariable('dradius','f4',('radius',),zlib=False)
            # drad.short_name        = 'DR'
            # drad.long_name   = 'aerosol size distribution bin sizes'
            # drad.units       = 'microns'
            # drad[:]          = self.DR*1e6

            # # Aerosol size distribution
            # for spc in ['TOT','BC','OC','DU','SS','SU']:
            #     varname = spc + 'dist'
            #     if spc == 'TOT':
            #         longname = 'total'
            #     else:
            #         longname = spc

            #     if hasattr(self,varname):
            #         dist = nc.createVariable(varname,'f4',('lev','radius',),zlib=zlib,fill_value=MISSING)
            #         dist.long_name     = '{} aerosol size distribution (dV/dlnr)'.format(longname)
            #         dist.missing_value = np.float32(MISSING)
            #         dist.units         = "microns^3/microns^2"
            #         dist[:]            = np.squeeze(self.__dict__[varname])

            # spc = 'colTOT'
            # varname = spc + 'dist'
            # dist = nc.createVariable(varname,'f4',('radius',),zlib=zlib,fill_value=MISSING)
            # dist.long_name     = 'column total aerosol size distribution (dV/dlnr)'
            # dist.missing_value = np.float32(MISSING)
            # dist.units         = "microns^3/microns^2"
            # dist[:]            = np.squeeze(self.__dict__[varname])


            # ang = nc.createVariable('angle','f4',('angle',),zlib=False)
            # ang.long_name      = 'scattering angle'
            # ang.units          = 'degrees'
            # ang[:]             = self.angle

            # p11 = nc.createVariable('P11','f4',('angle',))
            # p11[:] = self.p11

            # p12 = nc.createVariable('P12','f4',('angle',))
            # p12[:] = self.p12

            # p33 = nc.createVariable('P33','f4',('angle',))
            # p33[:] = self.p33   

            # p34 = nc.createVariable('P34','f4',('angle',))
            # p34[:] = self.p34                        

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

    outRoot = './graspConfig_12_Osku_DrySU_V1/'


    #albedoType   = 'BPDF'
    #albedoType = 'MODIS_BRDF_BPDF'
    #albedoType = 'LAMBERTIAN_BPDF'
    albedoType = 'GissCX'
    #albedoType = 'CX'
    #albedoType = 'LAMBERTIAN'
    #albedoType = 'MODIS_BRDF'

    aerosol        = True # true if you want aerosols in simulation
    simple_aerosol = True # true if exponential aerosol with constnat optical properties
    emde           = False # use pmom from Emde et al paper
    grasp          = False # fit legendre coefficient to GRASP phase function
    grasp_higher   = False # fit higher order spherical functions to GRASP phase function
    rayleigh       = False  # use rayleigh phase function
    dark           = True  #use if you want lambertian surface with albedo = 0
    nstreams       = 12
    plane_parallel = True


    channels  = 800, #865, #470,    # 410,440,470,550,670,865,1020,1650,2100  #
    for channel in channels:
        chd      = get_chd(channel)

        if (albedoType == 'LAMBERTIAN') and (dark is True):
            sfcname = 'nosurface'
        elif albedoType == 'MODIS_BRDF':
            sfcname = 'BRDF'
        else:
            sfcname = albedoType

        if aerosol is False:
            aname = 'rayleigh'
        else:
            if simple_aerosol:
                #aname = 'rayleigh+simple_aerosol'
                aname = 'rayleigh+simple_aerosol'
                if emde:
                    aname = 'rayleigh+simple_aerosol+emde'
                elif grasp:
                    aname = 'rayleigh+simple_aerosol+grasp'
                elif grasp_higher:
                    aname = 'simple_aerosol+grasp_higher'
            else:
                aname = 'rayleigh+aerosol'

        outDir    = '{}/benchmark_{}_{}'.format(outRoot,aname,sfcname)
        outFile   = '{}/calipso-g5nr.vlidort.vector.{}.{}.nc4'.format(outDir,albedoType,chd)
        
        rcFile   = 'rc/Aod_EOS.{}.rc'.format(channel)
        verbose  = True

        # Initialize VLIDORT class getting aerosol optical properties
        # -----------------------------------------------------------
        vlidort = BENCHMARK(albedoType,channel,outFile,rootDir,rcFile,
                                verbose=verbose,aerosol=aerosol,dark=dark,
                                simple_aerosol=simple_aerosol,emde=emde,
                                grasp=grasp,grasp_higher=grasp_higher,rayleigh=rayleigh,
                                nstreams=nstreams,plane_parallel=plane_parallel)

       
        vlidort.runVLIDORT()
