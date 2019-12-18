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
import VLIDORT_POLAR_ 
from copyvar  import _copyVar
from scipy import interpolate
from   MAPL  import config

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
SDS_MET = [] #[CLDTOT]
SDS_INV = ['FRLAND']
SDS_CX = ['U10M','V10M']

ncALIAS = {'LONGITUDE': 'trjLon',
           'LATITUDE': 'trjLat'}


nMom     = 300

VZAdic = {'POLDER': np.array([3.66, 11., 18.33, 25.66, 33, 40.33, 47.66, 55.0])}

HGTdic = {'LEO': 705,
          'ISS': 400}


SurfaceFuncs = {'MODIS_BRDF'     : 'readSampledMODISBRDF',
                'MODIS_BRDF_BPDF': 'readSampledMODISBRDF',
                'LAMBERTIAN'     : 'readSampledLER',
                'CX'             : 'readSampledWindCX'}

WrapperFuncs = {'MODIS_BRDF'     : VLIDORT_POLAR_.vector_brdf_modis,
                'MODIS_BRDF_BPDF': VLIDORT_POLAR_.vector_brdf_modis_bpdf,
                'LAMBERTIAN'     : VLIDORT_POLAR_.vector_lambert,
                'CX'             : VLIDORT_POLAR_.vector_gisscx}   

LandAlbedos  = 'MODIS_BRDF','MODIS_BRDF_BPDF','LAMBERTIAN','CX'


MISSING = -1.e+20

class POLAR_VLIDORT(object):
    """
    Everything needed for calling VLIDORT
    GEOS-5 has already been sampled on lidar track
    """
    def __init__(self,inFile,outFile,rcFile,albedoType,
                channel,VZA,hgtss,
                brdfFile=None,
                ndviFile=None,
                lcFile=None,
                lerFile=None,
                verbose=False,
                extOnly=False,
                distOnly=False,
                outFileDist=None):
        self.SDS_AER     = SDS_AER
        self.SDS_MET     = SDS_MET
        self.SDS_INV     = SDS_INV
        self.SDS_CX      = SDS_CX
        self.AERNAMES    = AERNAMES
        self.inFile      = inFile
        self.outFile     = outFile
        self.outFileDist = outFileDist
        self.albedoType  = albedoType
        self.rcFile      = rcFile
        self.channel     = channel
        self.verbose     = verbose
        self.nMom        = nMom
        self.brdfFile    = brdfFile
        self.lcFile      = lcFile
        self.ndviFile    = ndviFile
        self.lerFile     = lerFile

        if not extOnly:
            # initialize empty lists
            for sds in self.SDS_AER+self.SDS_MET+self.SDS_INV+self.SDS_CX:
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
            
            if not distOnly:
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

        self.iLand = self.FRLAND >= 0.99
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

        if len(self.SDS_MET) > 0:
            col = 'met_Nv'
            if self.verbose: 
                print 'opening file',self.inFile.replace('%col',col)        
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
            print 'opening file',self.inFile.replace('%col',col)        
        nc       = Dataset(self.inFile.replace('%col',col))

        for sds in self.SDS_CX:
            sds_ = sds
            if sds in ncALIAS:
                sds_ = ncALIAS[sds]
            var = nc.variables[sds_][:]
            self.__dict__[sds].append(var)

        for sds in self.SDS_CX:
            self.__dict__[sds] = np.concatenate(self.__dict__[sds])

 

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
        rmin      = 0.005e-6   #meters 
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
        rmax0     = 0.3e-6     #meters
        sigma     = 2.00
        rhop0     = 1000       # Density of dry particles [kg m-3]

        spc = 'BCPHOBIC'
        self.logNormalDistribution(spc,r0,rmin,rmax0,sigma,rhop0)
        spc = 'BCPHILIC'
        self.logNormalDistribution(spc,r0,rmin,rmax0,sigma,rhop0)

        # OC
        r0        = 0.0212e-6  #meters
        rmax0     = 0.3e-6     #meters
        sigma     = 2.20
        rhop0     = 1800       #[kg m-3]

        spc = 'OCPHOBIC'
        self.logNormalDistribution(spc,r0,rmin,rmax0,sigma,rhop0)
        spc = 'OCPHILIC'
        self.logNormalDistribution(spc,r0,rmin,rmax0,sigma,rhop0)

        # SU
        r0        = 0.0695e-6  #meters
        rmax0     = 0.3e-6     #meters
        sigma     = 2.03
        rhop0     = 1700       #[kg m-3]

        spc = 'SU'
        self.logNormalDistribution(spc,r0,rmin,rmax0,sigma,rhop0)

        # Dust
        # ----------
        self.dustDistribution()

        # Sea Salt
        # ------------
        self.seasaltDistribution()

        # Mixture Distribution
        # -------------
        self.TOTdist = self.BCPHILICdist + self.BCPHOBICdist + self.OCPHILICdist + self.OCPHOBICdist + self.SUdist + self.DUdist + self.SSdist
        self.BCdist  = self.BCPHILICdist + self.BCPHOBICdist
        self.OCdist  = self.OCPHILICdist + self.OCPHOBICdist

        # convert to dV/dlnR [microns^3/microns^2] 
        for i,r in enumerate(self.R):
            self.TOTdist[:,:,i] = self.TOTdist[:,:,i]*self.R[i]*1e6
            self.BCdist[:,:,i] = self.BCdist[:,:,i]*self.R[i]*1e6
            self.OCdist[:,:,i] = self.OCdist[:,:,i]*self.R[i]*1e6
            self.SUdist[:,:,i] = self.SUdist[:,:,i]*self.R[i]*1e6
            self.SSdist[:,:,i] = self.SSdist[:,:,i]*self.R[i]*1e6
            self.DUdist[:,:,i] = self.DUdist[:,:,i]*self.R[i]*1e6


    def getRMAX(self):
        cf = config.Config(self.rcFile)

        #lognormals
        spclist = 'BC','OC','SU'
        rmax0   = 0.3e-6
        RMAX    = 0
        for spc in spclist:
            # Read optics table            
            optable = cf('filename_optical_properties_{}'.format(spc))
            nc      = Dataset(optable)
            if spc == 'SU':
                gfTable = nc.variables['growth_factor'][0,:]
            else:
                gfTable = nc.variables['growth_factor'][1,:]

            rhTable = nc.variables['rh'][:]
            nc.close()
            fTable = interpolate.interp1d(rhTable, gfTable)            
            gf = fTable(rhTable[-1])
            rmax  = gf*rmax0
            if rmax > RMAX: RMAX = rmax

        # Sea-salt
        c1 = 0.7674
        c2 = 3.079
        c3 = 2.573e-11
        c4 = -1.424        
        rhUse = 0.95
        # read optics table
        optable = cf('filename_optical_properties_SS')
        nc      = Dataset(optable)
        rUp     = nc.variables['rUp'][:]
        nc.close()
        rMaxUse = rUp[-1]
        rMaxCM  = rMaxUse*100.
        rmax    = (c1*rMaxCM**c2 /(c3*rMaxCM**c4 - np.log10(rhUse))+rMaxCM**3.)**(1./3.)/100.
        if rmax > RMAX: RMAX = rmax

        # Dust
        # read optics table
        optable = cf('filename_optical_properties_DU')
        nc      = Dataset(optable)
        rUp     = nc.variables['rUp'][:]
        nc.close()
        rmax    = rUp[-1]
        if rmax > RMAX: RMAX = rmax

        return RMAX


    def seasaltDistribution(self):
        # Density of dry particles [kg m-3]
        rhop0 = 2200.

        #constants for adjusting size bins for RH
        c1 = 0.7674
        c2 = 3.079
        c3 = 2.573e-11
        c4 = -1.424        


        # master bins
        R    = self.R
        DR   = self.DR
        RLOW = self.RLOW
        RUP  = self.RUP

        # Read optics table
        cf = config.Config(self.rcFile)
        optable = cf('filename_optical_properties_SS')

        # Read in growth factors and rh tables
        nc = Dataset(optable)
        rhTable = nc.variables['rh'][:]
        gfTable = nc.variables['growth_factor'][:]

        # Major bins
        rMaxMaj = nc.variables['rUp'][:]
        rMinMaj = nc.variables['rLow'][:]

        nc.close()

        # number of major bins
        nbinMaj = len(rMaxMaj)

        #number of minor bins
        nbinMin = 40

        # set up empty column size distribution array
        nlev = 72
        SPCdist = np.empty([self.nobs,nlev,len(R)])

        # rh (relative humidity) dims are [ntyme,nlev]
        # mr (mixing ratio) dims are [ntyme,nlev]
        # air (air density) dims are [ntyme,nlev]
        rh = self.RH.copy()
        ss001,ss002,ss003,ss004,ss005 = self.SS001,self.SS002,self.SS003,self.SS004,self.SS005
        air = self.AIRDENS

        # put all sea-salt mixing ratios in one array for convenience
        SS = np.zeros([self.nobs,nlev,5])
        SS[:,:,0] = ss001
        SS[:,:,1] = ss002
        SS[:,:,2] = ss003
        SS[:,:,3] = ss004
        SS[:,:,4] = ss005


        # 0 <= rh <= 0.95
        # this is what is done in Chem_MieMod.F90
        rh[rh < 0] = 0
        rh[rh > 0.95] = 0.95   

        # loop throuhg time steps
        for t in range(self.nobs):    
            # loop through major bins
            for iBin in range(nbinMaj):
                # get the radii of the bin
                rmin = rMinMaj[iBin]
                rmax = rMaxMaj[iBin]
                rMinCM = rmin*100.
                rMaxCM = rmax*100.                    

                # get growth factors table for this bin
                fTable = interpolate.interp1d(rhTable, gfTable[iBin,:])
                gf = fTable(rh[t,:])  

                # loop through layers 
                # get size distribution for each layer
                for k in range(nlev):

                    #adjust bin edges for humidified particles
                    rhUse = rh[t,k]
                    rMinUse = (c1*rMinCM**c2 /(c3*rMinCM**c4 - np.log10(rhUse))+rMinCM**3.)**(1./3.)/100.                                      
                    rMaxUse = (c1*rMaxCM**c2 /(c3*rMaxCM**c4 - np.log10(rhUse))+rMaxCM**3.)**(1./3.)/100.  

                    # Determine the dNdr of the particle size distribution using the
                    # Gong 2003 particle sub-bin distribution
                    rrat = rmin/rMinUse
                    r80Rat = 1.65*rrat      # ratio of the r80 radius to the wet radius
                    r80  = R*r80Rat * 1.e6  # radius in r80 space in um
                    dr80 = DR*r80Rat * 1.e6

                    aFac = 4.7*(1.+30.*r80)**(-0.017*r80**(-1.44))
                    bFac = (0.433-np.log10(r80))/0.433
                    dndr80 = 1.373*r80**(-aFac)*(1.+0.057*r80**3.45)*10.**(1.607*np.exp(-bFac**2.))
                    dndr = dndr80 * r80Rat

                    # Truncate distribution according to rlow and rup
                    ii = RUP <= rMinUse
                    dndr[ii] = 0

                    ii = RLOW >= rMaxUse
                    dndr[ii] = 0

                    # deal with lowest bin
                    # number concentration is scaled to the 
                    # fraction of the bin that is covered by
                    # the aerosol distribution
                    bini = np.arange(len(R))
                    i = bini[RLOW < rMinUse][-1]
                    drtilda = RUP[i] - rMinUse
                    dndr[i] = dndr[i]*drtilda/DR[i]

                    #deal with the highest bin
                    # number concentration is scaled to the 
                    # fraction of the bin that is covered by
                    # the aerosol distribution                    
                    i = bini[RLOW < rMaxUse][-1]
                    drtilda = rMaxUse - RLOW[i]
                    dndr[i] = dndr[i]*drtilda/DR[i]

                    # Now get the volume distribution
                    # dvdr
                    dvdr = 4./3.*np.pi*R**3.*dndr

                    # Get aerosol DRY! volume concentration
                    mr = SS[t,k,iBin]
                    M0 = mr*air[t,k]
                    V0 = M0/rhop0

                    # Get the Wet volume
                    Vwet = V0*gf[k]**3

                    # normalize dvdr so the integral is equal to the 
                    # wet volume
                    dvdr = dvdr*Vwet/np.sum(dvdr*DR)

                    # add this aerosol distribution to the master
                    SPCdist[t,k,:] = SPCdist[t,k,:] + dvdr

        self.__dict__['SSdist'] = SPCdist
        

    def dustDistribution(self):
        # master bins
        R    = self.R
        DR   = self.DR
        RLOW = self.RLOW
        RUP  = self.RUP

        # read optics table
        cf = config.Config(self.rcFile)
        optable = cf('filename_optical_properties_DU')
        nc      = Dataset(optable)
        # Major bins
        rMaxMaj = nc.variables['rUp'][:]
        rMinMaj = nc.variables['rLow'][:]
        # density
        rhop0   = nc.variables['rhop'][:,0]

        nc.close()

        # number of major bins
        nbinMaj = len(rMaxMaj)

        #number of minor bins
        nbinMin = 10

        # set up empty column size distribution array
        nlev = 72
        SPCdist = np.zeros([self.nobs,nlev,len(R)])

        # mr (mixing ratio) dims are [ntyme,nlev]
        # air (air density) dims are [ntyme,nlev]
        du001,du002,du003,du004,du005  = self.DU001,self.DU002,self.DU003,self.DU004,self.DU005
        air = self.AIRDENS

        # put all dust mixing ratios in one array for convenience
        DU = np.empty([self.nobs,nlev,5])
        DU[:,:,0] = du001
        DU[:,:,1] = du002
        DU[:,:,2] = du003
        DU[:,:,3] = du004
        DU[:,:,4] = du005

        # loop through orbit
        for t in range(self.nobs):    
            # loop through layers
            for k in range(nlev):

                #do first bin
                iBin = 0
                nBinFirst = 4
                rMinFirst = np.array([0.1,0.18,0.3,0.6])*1.e-6
                rMaxFirst = np.array([0.18,0.3,0.6,1.])*1.e-6
                fMass = [0.009, 0.081, 0.234, 0.676] 

                for iBinfirst in range(nBinFirst):
                    # get the radii of the bin
                    rmin = rMinFirst[iBinfirst]
                    rmax = rMaxFirst[iBinfirst]

                    dndr = R**(-4.)

                    # Truncate distribution according to rlow and rup
                    ii = RUP <= rmin
                    dndr[ii] = 0

                    ii = RLOW >= rmax
                    dndr[ii] = 0

                    # deal with lowest bin
                    # number concentration is scaled to the 
                    # fraction of the bin that is covered by
                    # the aerosol distribution
                    bini = np.arange(len(R))
                    i = bini[RLOW < rmin][-1]
                    drtilda = RUP[i] - rmin
                    dndr[i] = dndr[i]*drtilda/DR[i]

                    #deal with the highest bin
                    # number concentration is scaled to the 
                    # fraction of the bin that is covered by
                    # the aerosol distribution                    
                    i = bini[RLOW < rmax][-1]
                    drtilda = rmax - RLOW[i]
                    dndr[i] = dndr[i]*drtilda/DR[i]

                    # Now get the volume distribution
                    # dvdr
                    # this is not quite right because
                    # dust is ellipsoid, but close enough
                    dvdr = 4./3.*np.pi*R**3.*dndr
                    
                    # Get aerosol DRY! volume concentration
                    mr = DU[t,k,iBin]
                    M0 = mr*air[t,k]
                    M0 = M0*fMass[iBinfirst]
                    V0 = M0/rhop0[iBin]

                    # Get the Wet volume
                    # same as dry because dust is hydrophobic
                    Vwet = V0

                    # normalize dvdr so the integral is equal to the 
                    # wet volume
                    dvdr = dvdr*Vwet/np.sum(dvdr*DR)

                    # add this aerosol distribution to the master
                    SPCdist[t,k,:] = SPCdist[t,k,:] + dvdr



                # do last 4 bins - they are easy
                for iBin in range(1,nbinMaj):
                    # get the radii of the bin
                    rmin = rMinMaj[iBin]
                    rmax = rMaxMaj[iBin]

                    dndr = R**(-4.)

                    # Truncate distribution according to rlow and rup
                    ii = RUP <= rmin
                    dndr[ii] = 0

                    ii = RLOW >= rmax
                    dndr[ii] = 0

                    # deal with lowest bin
                    # number concentration is scaled to the 
                    # fraction of the bin that is covered by
                    # the aerosol distribution
                    bini = np.arange(len(R))
                    i = bini[RLOW < rmin][-1]
                    drtilda = RUP[i] - rmin
                    dndr[i] = dndr[i]*drtilda/DR[i]

                    #deal with the highest bin
                    # number concentration is scaled to the 
                    # fraction of the bin that is covered by
                    # the aerosol distribution                    
                    i = bini[RLOW < rmax][-1]
                    drtilda = rmax - RLOW[i]
                    dndr[i] = dndr[i]*drtilda/DR[i]
        

                    # Now get the volume distribution
                    # dvdr
                    # this is not quite right because
                    # dust is ellipsoid, but close enough
                    dvdr = 4./3.*np.pi*R**3.*dndr
                    
                    # Get aerosol DRY! volume concentration
                    mr = DU[t,k,iBin]
                    M0 = mr*air[t,k]
                    V0 = M0/rhop0[iBin]

                    # Get the Wet volume
                    # same as dry because dust is hydrophobic
                    Vwet = V0

                    # normalize dvdr so the integral is equal to the 
                    # wet volume
                    dvdr = dvdr*Vwet/np.sum(dvdr*DR)

                    # add this aerosol distribution to the master
                    SPCdist[t,k,:] = SPCdist[t,k,:] + dvdr

        self.__dict__['DUdist'] = SPCdist



    def logNormalDistribution(self,spc,r0,rmin,rmax0,sigma,rhop0):

        # master bins
        R    = self.R
        DR   = self.DR
        RLOW = self.RLOW
        RUP  = self.RUP

        # Read optics table
        cf = config.Config(self.rcFile)
        optable = cf('filename_optical_properties_{}'.format(spc[0:2]))

        # Read in growth factors and rh tables
        nc = Dataset(optable)
        rhTable = nc.variables['rh'][:]
        if 'PHOBIC' in spc:
            gfTable = nc.variables['growth_factor'][0,:]
        elif 'PHILIC' in spc:
            gfTable = nc.variables['growth_factor'][1,:]
        else:
            gfTable = nc.variables['growth_factor'][0,:]

        nc.close()
        fTable = interpolate.interp1d(rhTable, gfTable)

        # set up empty column size distribution array
        nlev = 72
        SPCdist = np.empty([self.nobs,nlev,len(R)])
        if spc == 'SU': 
            spc_ = 'SO4'
        else:
            spc_ = spc

        # rh dims are [ntyme,nlev]
        # mr (mixing ratio) dims are [ntyme,nlev]
        # air (air density) dims are [ntyme,nlev]
        rh  = self.RH.copy()
        mr  = self.__dict__[spc_]
        air = self.AIRDENS

        # 0 <= rh <= 0.99
        # this is what is done in Chem_MieMod.F90
        rh[rh < 0] = 0
        rh[rh > 0.99] = 0.99   

        # loop throuhg orbit
        for t in range(self.nobs):        
            gf = fTable(rh[t,:])

            #  Now create the particle properties for this humidified particle
            rmode = gf*r0
            rmax  = gf*rmax0   

            # loop through layers 
            # get size distribution for each layer
            for k in range(nlev):
                # Get the bins
                # r, dr, rlow, rup = self.logBins(rmin,rmax[k])

                # Now get the aeorosl number distribution
                # dndr
                rNum = rmode[k]
                lsigma = np.log(sigma)
                C      = np.sqrt(2.*np.pi)
                dndr = (1./(R*lsigma*C))*np.exp(-(np.log(R/rNum)**2.)/(2.*lsigma**2.)) 

                # Truncate distribution according to rlow and rup
                ii = RUP <= rmin
                dndr[ii] = 0

                ii = RLOW >= rmax[k]
                dndr[ii] = 0

                # deal with lowest bin
                # number concentration is scaled to the 
                # fraction of the bin that is covered by
                # the aerosol distribution
                bini = np.arange(len(R))
                i = bini[RLOW < rmin][-1]
                drtilda = RUP[i] - rmin
                dndr[i] = dndr[i]*drtilda/DR[i]

                #deal with the highest bin
                # number concentration is scaled to the 
                # fraction of the bin that is covered by
                # the aerosol distribution                    
                i = bini[RLOW < rmax[k]][-1]
                drtilda = rmax[k] - RLOW[i]
                dndr[i] = dndr[i]*drtilda/DR[i]

                # Now get the volume distribution
                # dvdr
                dvdr = 4./3.*np.pi*R**3.*dndr

                
                # Get aerosol DRY! volume concentration
                M0 = mr[t,k]*air[t,k]
                V0 = M0/rhop0

                # Get the Wet volume
                Vwet = V0*gf[k]**3

                # normalize dvdr so the integral is equal to the 
                # wet volume                    
                dvdr = dvdr*Vwet/np.sum(dvdr*DR)

                # # use the fact that for a log-normal distribution
                # # volume = 4/3*pi*N0*rmode^3*exp(9/2*ln(sigma)^2)
                # C = (4./3.)*np.pi*np.exp((9./2.)*lsigma**2)
                # N0   = Vwet/(C*rNum**3)

                # add this aerosol distribution to the master
                SPCdist[t,k,:] = dvdr

        self.__dict__[spc+'dist'] = SPCdist


    def logBins(self,rmin,rmax,nbin=None):
        """ Get size bins for lognormal distribution """
        # From P. Colarco's runmie_lognormal.pro
        # Now we want to set up the aerosol bins.  
        # We impose a cutoff on the maximum radius to consider.
        # We consider a minimum radius of 0.005 um (following GADS) 
        # and use a number of bins equal to at
        # least 20 bins per decade of the ratio rmax/rmin.

        # The bins are centered in volume betwen rlow and rup
        # That is...r^3-rlow^3 = rup^3-r^3
        # rup^3 = rlow^3*rmrat, which can be solved to find r given rmrat
        # and a desired rlow, e.g. rmin = 1.d-6*((1.+rmrat)/2.)^(1.d/3)
        # where the desired rlow = 1.d-6 in this example.
        # The meaning of r is that it is the radius of the particle with 
        # the average volume of the bin.        

        if nbin is None:
            decade = np.floor(np.log10(rmax/rmin)+1)
            nbin   = int(20*decade)
        rat    = rmax/rmin
        rmRat  = (rat**3.0)**(1.0/nbin)
        rMinUse = rmin*((1.+rmRat)/2.)**(1.0/3.0)  

        cpi = np.pi*4./3.
        rvolmin = cpi*rMinUse**3.
        vrfact = ( (3./2./np.pi / (rmRat+1))**(1./3.))*(rmRat**(1./3.) - 1.)  

        rvol     = np.zeros(nbin)
        rvolup   = np.zeros(nbin)
        r        = np.zeros(nbin)
        rup      = np.zeros(nbin)
        dr       = np.zeros(nbin)
        rlow     = np.zeros(nbin)
        for ibin in range(nbin):
            rvol[ibin]   = rvolmin*rmRat**float(ibin)
            rvolup[ibin] = 2.*rmRat/(rmRat+1.)*rvol[ibin]
            r[ibin]       = (rvol[ibin]/cpi)**(1./3.)
            rup[ibin]     = (rvolup[ibin]/cpi)**(1./3.)
            dr[ibin]      = vrfact*rvol[ibin]**(1./3.)
            rlow[ibin]    = rup[ibin] - dr[ibin]
        
        return r, dr, rlow, rup


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
            print 'opening BRDF file ',self.brdfFile
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
            print 'opening BRDF abledo file ',self.brdfFile
        nc = Dataset(self.brdfFile)

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
                      " --intensive"  +\
                      " --du --su --ss --oc --bc"    
                      

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
        surfList = []
        if self.nobsLand > 0:
            self.iLand = np.arange(len(self.iGood))[self.iGood & self.iLand]
            surfList.append('Land')
        if self.nobsSea > 0:
            self.iSea = np.arange(len(self.iGood))[self.iGood & self.iSea]
            surfList.append('Sea')
        self.iGood = np.arange(len(self.iGood))[self.iGood]

        # Initiate output arrays
        nangles = len(self.VZA)
        ntime   = len(self.tyme)
        nlev    = self.tau.shape[0]
        self.I = np.ones([ntime,2*nangles])*MISSING
        self.Q = np.ones([ntime,2*nangles])*MISSING
        self.U = np.ones([ntime,2*nangles])*MISSING
        self.reflectance = np.ones([ntime,2*nangles])*MISSING
        self.surf_reflectance = np.ones([ntime,2*nangles])*MISSING
        self.BR_Q = np.ones([ntime,2*nangles])*MISSING
        self.BR_U = np.ones([ntime,2*nangles])*MISSING
        self.ROT = np.ones([ntime,nlev])*MISSING

        #backward directions, forward directions
        VZA = np.append(self.VZA[::-1],self.VZA)

        # loop though LAND and SEA
        for surface in surfList:

            iGood = self.__dict__['i'+surface]
            sza  = self.SZA[iGood]
            raaf = self.RAAf[iGood]
            raab = self.RAAb[iGood]

            raaf.shape = raaf.shape + (1,)
            raaf = np.repeat(raaf,len(self.VZA),axis=1)
            
            raab.shape = raab.shape + (1,)        
            raab = np.repeat(raab,len(self.VZA),axis=1)

            RAA = np.append(raab,raaf,axis=1)
            RAA[RAA < 0] = RAA[RAA<0]+360.0

            tau = self.tau[:,:,iGood]
            ssa = self.ssa[:,:,iGood]
            pmom = self.pmom[:,:,iGood,:,:]
            pe   = self.pe[:,iGood]
            ze   = self.ze[:,iGood]
            te   = self.te[:,iGood]

            if surface == 'Land':        
                albedoType = self.albedoType

            else:
                albedoType = 'CX'

            # Get VLIDORT wrapper name from dictionary
            vlidortWrapper = WrapperFuncs[albedoType]

            # loop through viewing angles and run VLIDORT
            for i in range(2*nangles):
                vza  = np.array([VZA[i]]*self.nobs)
                raa  = RAA[:,i]

                runFunction = self.__getattribute__(albedoType+'_run')
                I,Q,U,reflectance,surf_reflectance,BR_Q,BR_U,ROT = runFunction(vlidortWrapper,tau,ssa,pmom,pe,ze,te,sza,raa,vza,iGood)
                                                    
                if i == 0:
                    self.ROT[iGood,:] = np.squeeze(ROT).T
                self.I[iGood,i] = np.squeeze(I)
                self.reflectance[iGood,i] = np.squeeze(reflectance)
                self.surf_reflectance[iGood,i] = np.squeeze(surf_reflectance)
                self.Q[iGood,i] = np.squeeze(Q)
                self.U[iGood,i] = np.squeeze(U) 
                self.BR_Q[iGood,i] = np.squeeze(BR_Q)
                self.BR_U[iGood,i] = np.squeeze(BR_U) 


        self.writeNC()
    #---
    def MODIS_BRDF_run(self,vlidortWrapper,tau,ssa,pmom,pe,ze,te,sza,raa,vza,iGood):
        kernel_wt = self.kernel_wt[:,:,iGood]
        param     = self.RTLSparam[:,:,iGood]                
        

        args = [self.channel,tau, ssa, pmom, 
                pe, ze, te, 
                kernel_wt, param, 
                sza, raa, vza, 
                MISSING,
                self.verbose]

        # Call VLIDORT wrapper function
        I, reflectance, ROT, surf_reflectance, Q, U, BR_Q, BR_U, rc = vlidortWrapper(*args)                        

        return I,Q,U,reflectance,surf_reflectance,BR_Q,BR_U,ROT        
    #---    
    def MODIS_BRDF_BPDF_run(self,vlidortWrapper,tau,ssa,pmom,pe,ze,te,sza,raa,vza,iGood):
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
        ROT              = np.zeros([nlev,self.nobsLand,1])

        for p in range(self.nobsLand):
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
                        sza[p:p+1], 
                        raa[p:p+1], 
                        vza[p:p+1], 
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

        return I,Q,U,reflectance,surf_reflectance,BR_Q,BR_U,ROT
    #---
    def LAMBERTIAN_run(self,vlidortWrapper,tau,ssa,pmom,pe,ze,te,sza,raa,vza,iGood):
        albedo = self.albedo[iGood,:]
        
        args = [self.channel,tau, ssa, pmom, 
                pe, ze, te, 
                albedo, 
                sza, raa, vza, 
                MISSING,
                self.verbose]

        # Call VLIDORT wrapper function
        I, reflectance, ROT, Q, U, rc = vlidortWrapper(*args)  
        surf_reflectance = albedo            

        BR_Q = None
        BR_U = None
        return I,Q,U,reflectance,surf_reflectance,BR_Q,BR_U,ROT
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

        sref = nc.createVariable('surf_reflectance_Q','f4',('time','view_angles',),zlib=zlib,fill_value=MISSING)
        sref.standard_name = '%.2f nm Surface Reflectance Q' %self.channel
        sref.long_name     = '%.2f nm Bi-Directional Surface Reflectance Q' %self.channel
        sref.missing_value = MISSING
        sref.units         = "None"
        sref[:]            = self.BR_Q

        sref = nc.createVariable('surf_reflectance_U','f4',('time','view_angles',),zlib=zlib,fill_value=MISSING)
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

        # Close the file
        # --------------
        nc.close()

        if self.verbose:
            print " <> wrote %s"%(self.outFile)

    #---
    def writeNCdist (self,zlib=True):
        """
        Write a NetCDF file of aerosol size distribution
        """
        km = 72

        if not os.path.exists(os.path.dirname(self.outFileDist)):
            os.makedirs(os.path.dirname(self.outFileDist))

        # Open NC file
        # ------------
        nc = Dataset(self.outFileDist,'w',format='NETCDF4_CLASSIC')

        # Set global attributes
        # ---------------------
        nc.title = 'GEOS-5 - Aerosol size distribution'
        nc.institution = 'NASA/Goddard Space Flight Center'
        nc.source = 'Global Model and Assimilation Office'
        nc.history = 'Aerosol size distribution of LEO satellite sampled on GEOS-5'
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
        nr = nc.createDimension('radius',len(self.R))

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

        rad = nc.createVariable('radius','f4',('radius',),zlib=False)
        rad.long_name   = 'aerosol radius'
        rad.units       = 'microns'
        rad[:]          = self.R*1e6

        # Aerosol size distribution
        for spc in ['TOT','BC','OC','DU','SS','SU']:
            varname = spc + 'dist'
            if spc == 'TOT':
                longname = 'total'
            else:
                longname = spc
            dist = nc.createVariable(varname,'f4',('time','lev','radius',),zlib=zlib,fill_value=MISSING)
            dist.long_name     = '{} aerosol size distribution (dV/dlnr)'.format(longname)
            dist.missing_value = MISSING
            dist.units         = "microns^3/microns^2"
            dist[:]            = self.__dict__[varname]

        # Close the file
        # --------------
        nc.close()

        if self.verbose:
            print " <> wrote %s"%(self.outFileDist)
    

def get_chd(channel):
    chd = '%.2f'%channel
    chd = chd.replace('.','d')

    return chd
    
#------------------------------------ M A I N ------------------------------------

if __name__ == "__main__":
    date     = datetime(2006,8,01,00)
    nymd     = str(date.date()).replace('-','')
    hour     = str(date.hour).zfill(2)
    format   = 'NETCDF4_CLASSIC'

    #rootDir  = '/discover/nobackup/projects/gmao/osse2/pub/c1440_NR/OBS/POLAR_LIDAR/CALIPSO/'
    rootDir  = '/nobackup/3/pcastell/POLAR_LIDAR/CALIPSO/'

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

    channel   = 470
    chd       = get_chd(channel)
    outDir    = '{}/LevelC/Y{}/M{}'.format(rootDir,date.year,str(date.month).zfill(2))
    outFile   = '{}/calipso-g5nr.vlidort.vector.{}.{}.{}_{}z_{}nm.nc4'.format(outDir,albedoType,VZAname,nymd,hour,chd)
    
    rcFile   = 'Aod_EOS.rc'
    orbit    = 'LEO'
    verbose  = True

    # Initialize VLIDORT class getting aerosol optical properties
    # -----------------------------------------------------------
    vlidort = POLAR_VLIDORT(inFile,outFile,rcFile,
                            albedoType,
                            channel,
                            VZAdic[VZAname],
                            HGTdic[orbit],
                            brdfFile=brdfFile,
                            ndviFile=ndviFile,
                            lcFile=lcFile,
                            verbose=verbose,
                            distOnly=False)

   
    # Run ext_sampler
    # vlidort.runExt()

    # Run VLIDORT
    # if vlidort.nobs > 0:
    #     vlidort.runVLIDORT()
