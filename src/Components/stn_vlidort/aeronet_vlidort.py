#!/usr/bin/env python

"""
    Calculates polarized BOA radiance at ground stations.
    Model fields have already been sampled using stn_sampler

    Adapted from polar_vlidort.py
    Patricia Castellanos, June, 2017

"""

import os
import MieObs_
from   netCDF4         import Dataset
from   mieobs          import  getAOPvector, getEdgeVars, getAOPint
import numpy           as np

from datetime          import datetime, timedelta
from dateutil.parser   import parse         as isoparser

from MAPL.constants    import *
import stnAngles_    
import VLIDORT_STN_ 
from scipy             import interpolate
import multiprocessing
from   stn_vlidort_aux import _copyVar, extrap1d, MieVARS, get_chd 
from   MAPL  import config

# Generic Lists of Varnames and Units
VNAMES_DU = ['DU001','DU002','DU003','DU004','DU005']
VNAMES_SS = ['SS001','SS002','SS003','SS004','SS005']
VNAMES_BC = ['BCPHOBIC','BCPHILIC']
VNAMES_OC = ['OCPHOBIC','OCPHILIC']
VNAMES_SU = ['SO4']
MieVarsNames = ['ext','scatext','backscat','aback_sfc','aback_toa','depol','ext2back','tau','ssa','g']
MieVarsUnits = ['km-1','km-1','km-1 sr-1','sr-1','sr-1','unitless','sr','unitless','unitless','unitless']

META    = ['DELP','PS','RH','AIRDENS','LONGITUDE','LATITUDE','isotime','stnName']
AERNAMES = VNAMES_SU + VNAMES_SS + VNAMES_OC + VNAMES_BC + VNAMES_DU
SDS_AER = META + AERNAMES
SDS_MET = []#['CLDTOT']
SDS_INV = ['FRLAND']

ncALIAS = {'LONGITUDE': 'stnLon',
           'LATITUDE' : 'stnLat'}

nMom     = 300

SurfaceFuncs = {'MODIS_BRDF'     : 'readSampledMODISBRDF',
                'MODIS_BRDF_BPDF': 'readSampledMODISBRDF',
                'LAMBERTIAN'     : 'readSampledLER'}

WrapperFuncs = {'MODIS_BRDF'     : VLIDORT_STN_.vector_brdf_modis,
                'MODIS_BRDF_BPDF': VLIDORT_STN_.vector_brdf_modis_bpdf,
                'LAMBERTIAN'     : VLIDORT_STN_.vector_lambert}   

LandAlbedos  = 'MODIS_BRDF','MODIS_BRDF_BPDF','LAMBERTIAN'


MISSING = -1.e+20

pp_angles = [  -6,  -5,  -4, -3.5,   -3, -2.5,  -2,  
                0,   2, 2.5,    3,  3.5,    4,   5,  
                6,   8,  10,   12,   14,   16,  20,   
               25,  30,  35,   40,   45,   50,  55, 
               60,  65,  70,   80,   90,  100, 110,  
              120, 130, 140,  150]

al_angles = [   0,  -6,   -5,   -4,  -3.5,    -3,  -2.5,
               -2,   2,  2.5,    3,   3.5,     4,     5,
                6,   7,    8,   10,    12,    14,    16,
               18,  20,   25,   30,    35,    40,    45,
               50,  60,   70,   80,    90,   100,   120,
              140, 160,  180, -180,  -160,  -140,  -120,
             -100, -90,  -80,  -70,   -60,   -50,   -45,
              -40, -35,  -30,  -25,   -20,   -18,   -16,
              -14, -12,  -10,   -8,    -7,    -6,    -6,
               -5,  -4, -3.5,   -3,  -2.5,    -2,     2,
              2.5,   3,  3.5,    4,     5,     6]              


pp_angles = [   0,   2, 2.5,    3,  3.5,    4,   5,  
                6,   8,  10,   12,   14,   16,  20,   
               25,  30,  35,   40,   45,   50,  55, 
               60,  65,  70,   80,   90,  100, 110,  
              120, 130, 140,  150]

al_angles = [   0,   2,  2.5,    3,   3.5,     4,     5,
                6,   7,    8,   10,    12,    14,    16,
               18,  20,   25,   30,    35,    40,    45,
               50,  60,   70,   80,    90,   100,   120,
              140, 160,  180]              


aeronet_r = [0.05,        0.065604,    0.086077,    0.112939,    0.148184,
             0.194429,    0.255105,    0.334716,    0.439173,    0.576227,
             0.756052,    0.991996,    1.301571,    1.707757,    2.240702,    
             2.939966,    3.857452,    5.06126,     6.640745,    8.713145,    
             11.432287,   15]  # microns

aeronet_r = np.array(aeronet_r)*1e-6  # meters

def unwrap_self_doMie(arg, **kwarg):
    return AERONET_VLIDORT.doMie(*arg, **kwarg)
 

class AERONET_VLIDORT(object):
    """
    Everything needed for calling VLIDORT
    GEOS-5 has already been sampled on lidar track
    """
    def __init__(self,inFile,rcFile,channel,
                invFile=None,
                outFileal=None,
                outFilepp=None,
                outFileext=None,
                outFileadd=None,
                albedoType=None,
                brdfFile=None,
                ndviFile=None,
                lcFile=None,
                lerFile=None,
                verbose=False,
                extOnly=False,
                extcol='aer_Nv',
                aeronet_r=aeronet_r):
        self.SDS_AER = SDS_AER
        self.SDS_MET = SDS_MET
        self.SDS_INV = SDS_INV
        self.AERNAMES = AERNAMES
        self.inFile  = inFile
        self.invFile = invFile
        self.outFileal = outFileal
        self.outFilepp = outFilepp
        self.outFileext = outFileext
        self.outFileadd = outFileadd
        self.albedoType = albedoType
        self.rcFile  = rcFile
        self.channel = channel
        self.verbose = verbose
        self.nMom    = nMom
        self.brdfFile = brdfFile        
        self.lcFile = lcFile
        self.ndviFile = ndviFile
        self.lerFile  = lerFile
        self.al_angles = al_angles
        self.pp_angles = pp_angles
        self.nal       = len(self.al_angles)
        self.npp       = len(self.pp_angles)
        self.aeronet_r = aeronet_r
        self.extcol    = extcol

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
        self.nstations = len(self.LONGITUDE)
        self.ntyme    = len(self.tyme)
        self.nobs  = [self.nstations]*self.ntyme
        self.iGood = [np.ones([self.nstations]).astype(bool)]*self.ntyme

        if not extOnly:
            if self.outFileadd is None:
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

            if self.outFileadd is not None:
                # Calculate column aerosol size distributions
                self.sizeDistribution()

                # get 440-870 angstrom exponent
                self.computeAE()

            # Calculate Scene Geometry
            self.calcAngles()

            if any(self.nobs) and (self.outFileadd is None):
                #Land-Sea Mask
                self.LandSeaMask()


        self.nobs = np.array(self.nobs)

    # --
    def BPDFinputs(self):
        """
        Read in NDVI and Landuse Coefficient- should have been sampled already
        """
        if self.verbose:
            print 'opening file',self.ndviFile

        nc = Dataset(self.ndviFile)
        NDVI = np.array(nc.variables['NDVI'][:])
        missing_value = nc.variables['NDVI'].missing_value
        I = NDVI < -900
        NDVI[I] = MISSING
        I = NDVI == missing_value
        NDVI[I] = MISSING
        nc.close()

        if self.verbose:
            print 'opening file',self.lcFile
        nc = Dataset(self.lcFile)
        BPDFcoef = nc.variables['BPDFcoef'][:]
        I = BPDFcoef < -900
        BPDFcoef[I] = MISSING
        nc.close()

        for s in range(self.ntyme):
            self.iGood[s] = self.iGood[s] & (NDVI[:,s] != MISSING)
            self.nobs[s]  = np.sum(self.iGood[s])

        #BPDFparam(nparam,nch,nobs)
        self.BPDFparam = [np.zeros([3,1,self.nstations])]*self.ntyme
        for t in range(self.ntyme):
            self.BPDFparam[t][0,0,:] = 1.5
            self.BPDFparam[t][1,0,:] = NDVI[:,t]
            self.BPDFparam[t][2,0,:] = BPDFcoef[:,t]        

    # ---
    def LandSeaMask(self):
        """
        Read in invariant dataset
        """
        col = 'asm_Nx'
        if self.verbose: 
            print 'opening file',self.invFile.replace('%col',col)
        nc       = Dataset(self.invFile.replace('%col',col))

        for sds in self.SDS_INV:
            sds_ = sds
            if sds in ncALIAS:
                sds_ = ncALIAS[sds]
            var = nc.variables[sds_][:]
            self.__dict__[sds].append(var)   

        # Make lists into arrays
        for sds in self.SDS_INV:
            self.__dict__[sds] = np.concatenate(self.__dict__[sds])            

        for t in range(self.ntyme):
            if self.albedoType in LandAlbedos:
                iGood = self.FRLAND[:,0] >= 0.99
            else:
                iGood = self.FRLAND[:,0] < 0.99

            self.iGood[t] = self.iGood[t] & iGood
            self.nobs[t]  = np.sum(self.iGood[t])

    # --
    def computeAtmos(self):
        self.pe = []
        self.ze = []
        self.te = []
        NAMES = ['AIRDENS','DELP','PS']
        for t in range(self.ntyme):
            inVars = MieVARS()
            for sds in NAMES:
                var = self.__dict__[sds]
                if len(var.shape) == 3:
                    inVars.__dict__[sds] = var[:,t,:]
                elif len(var.shape) == 2:
                    inVars.__dict__[sds] = var[:,t]

            pe, ze, te = getEdgeVars(inVars)

            self.pe.append(pe) # (km,nobs)
            self.ze.append(ze) 
            self.te.append(te)

    # --
    def calcAngles(self):
        self.SZA   = []
        self.SAA   = []
        for t in self.tyme:
            SZA = []
            SAA = []
            for s in range(self.nstations):
                CLAT = self.LATITUDE[s]
                CLON = self.LONGITUDE[s]
                
                year  = t.year 
                month = t.month 
                day   = t.day 
                hour  = t.hour 
                minute = t.minute 
                second = t.second             

                saa, sza = stnAngles_.sunangles(year,month,day,hour,minute,second,
                                                  CLAT,CLON,
                                                  0.0)
                SZA.append(sza[0])
                SAA.append(saa[0])

            self.SZA.append(np.array(SZA))
            self.SAA.append(np.array(SAA))

        # Limit SZAs
        for t in range(self.ntyme):
            iGood = self.SZA[t] < 80
            self.iGood[t] = self.iGood[t] & iGood
            self.nobs[t] = np.sum(self.iGood[t])     

            iGood = self.SZA[t] >= 80
            self.SZA[t][iGood] = MISSING
            self.SAA[t][iGood] = MISSING

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


        # # try higher res bins
        # dr     = 0.01 
        # RMIN   = RLOW[0]
        # RMAX   = RUP[-1]
        # deltaR = np.log(RMAX) - np.log(RMIN)
        # nbin   = int(np.ceil(deltaR/dr))
        # R = np.zeros(nbin)
        # DR = np.zeros(nbin)
        # RLOW = np.zeros(nbin)
        # RUP  = np.zeros(nbin)

        # RLOW[0] = RMIN
        # RUP[0]  = np.exp(np.log(RLOW[0]) + dr)
        # R[0]    = np.exp(np.log(RLOW[0]) + 0.5*dr)
        # for i in range(1,nbin):
        #     RLOW[i] = np.exp(np.log(RLOW[i-1]) + dr)
        #     RUP[i]  = np.exp(np.log(RLOW[i]) + dr)
        #     R[i]    = np.exp(np.log(RLOW[i]) + 0.5*dr)

        # self.R = R
        # self.DR = RUP - RLOW
        # self.RLOW = RLOW
        # self.RUP  = RUP



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
        self.TOTdist = []
        for st in range(self.nstations):
            TOTdist = self.BCPHILICdist[st] 
            TOTdist = TOTdist + self.BCPHOBICdist[st] 
            TOTdist = TOTdist + self.OCPHILICdist[st] 
            TOTdist = TOTdist + self.OCPHOBICdist[st] 
            TOTdist = TOTdist + self.SUdist[st] 
            TOTdist = TOTdist + self.DUdist[st] 
            TOTdist = TOTdist + self.SSdist[st]

            self.TOTdist.append(TOTdist)

        # Interpolate to AERONET radius
        self.AERdist = []
        for st in range(self.nstations):
            TOTdist = self.TOTdist[st]
            AERdist = np.zeros([self.ntyme,len(self.aeronet_r)])
            for t in range(self.ntyme):
                fTable = interpolate.interp1d(self.R,TOTdist[t,:])
                AERdist[t,:]   = fTable(self.aeronet_r)*self.aeronet_r*1e6

            self.AERdist.append(AERdist)


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

        # initiate output variables
        self.__dict__['SSdist'] = []

        # loop through stations
        # mr (mixing ratio) dims are [ntyme,nlev]  
        nlev = 72      
        for ss001,ss002,ss003,ss004,ss005,rh,air in zip(self.SS001,self.SS002,self.SS003,self.SS004,self.SS005,self.RH,self.AIRDENS):
            # set up empty column size distribution array
            SPCdist = np.zeros((self.ntyme,)+R.shape)

            # put all sea-salt mixing ratios in one array for convenience
            SS = np.empty([self.ntyme,nlev,5])
            SS[:,:,0] = ss001
            SS[:,:,1] = ss002
            SS[:,:,2] = ss003
            SS[:,:,3] = ss004
            SS[:,:,4] = ss005


            # 0 <= rh <= 0.95
            # this is what is done in Chem_MieMod.F90
            rh = rh.copy()
            rh[rh < 0] = 0
            rh[rh > 0.95] = 0.95   

            # loop throuhg time steps
            for t in range(self.ntyme):    
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

                        # get the bins
                        r, dr, rlow, rup = self.logBins(rMinUse,rMaxUse,nbin=nbinMin)

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
                        SPCdist[t,:] = SPCdist[t,:] + dvdr

            self.__dict__['SSdist'].append(SPCdist)
        

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

        # initiate output variables
        self.__dict__['DUdist'] = []

        # loop through stations
        # mr (mixing ratio) dims are [ntyme,nlev]  
        nlev = 72      
        for du001,du002,du003,du004,du005,air  in zip(self.DU001,self.DU002,self.DU003,self.DU004,self.DU005,self.AIRDENS):
            # set up empty column size distribution array
            SPCdist = np.zeros((self.ntyme,)+R.shape)

            # put all dust mixing ratios in one array for convenience
            DU = np.empty([self.ntyme,nlev,5])
            DU[:,:,0] = du001
            DU[:,:,1] = du002
            DU[:,:,2] = du003
            DU[:,:,3] = du004
            DU[:,:,4] = du005

            # loop through time steps
            for t in range(self.ntyme):    

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
                    r, dr, rlow, rup = self.logBins(rmin,rmax,nbin=nbinMin)

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
                    mr = DU[t,:,iBin]
                    M0 = mr*air[t,:]
                    M0 = M0*fMass[iBinfirst]
                    V0 = M0.sum()/rhop0[iBin]

                    # Get the Wet volume
                    # same as dry because dust is hydrophobic
                    Vwet = V0

                    # normalize dvdr so the integral is equal to the 
                    # wet volume
                    dvdr = dvdr*Vwet/np.sum(dvdr*DR)

                    # # use the fact that for a log-normal distribution
                    # # volume = 4/3*pi*N0*rmode^3*exp(9/2*ln(sigma)^2)
                    # C = (4./3.)*np.pi*np.exp((9./2.)*lsigma**2)
                    # N0   = Vwet/(C*rNum**3)

                    # add this aerosol distribution to the master
                    SPCdist[t,:] = SPCdist[t,:] + dvdr



                # do last 4 bins - they are easy
                for iBin in range(1,nbinMaj):
                    # get the radii of the bin
                    rmin = rMinMaj[iBin]
                    rmax = rMaxMaj[iBin]
                    r, dr, rlow, rup = self.logBins(rmin,rmax,nbin=nbinMin)

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
                    mr = DU[t,:,iBin]
                    M0 = mr*air[t,:]
                    V0 = M0.sum()/rhop0[iBin]

                    # Get the Wet volume
                    # same as dry because dust is hydrophobic
                    Vwet = V0

                    # normalize dvdr so the integral is equal to the 
                    # wet volume
                    dvdr = dvdr*Vwet/np.sum(dvdr*DR)

                    # # use the fact that for a log-normal distribution
                    # # volume = 4/3*pi*N0*rmode^3*exp(9/2*ln(sigma)^2)
                    # C = (4./3.)*np.pi*np.exp((9./2.)*lsigma**2)
                    # N0   = Vwet/(C*rNum**3)

                    # add this aerosol distribution to the master
                    SPCdist[t,:] = SPCdist[t,:] + dvdr

            self.__dict__['DUdist'].append(SPCdist)



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

        # loop through stations
        # rh dims are [ntyme,nlev]
        # mr (mixing ratio) dims are [ntyme,nlev]
        nlev = 72
        self.__dict__[spc+'dist'] = []
        if spc == 'SU': 
            spc_ = 'SO4'
        else:
            spc_ = spc

        for rh,mr,air  in zip(self.RH,self.__dict__[spc_],self.AIRDENS):
            # 0 <= rh <= 0.99
            # this is what is done in Chem_MieMod.F90
            rh = rh.copy()
            rh[rh < 0] = 0
            rh[rh > 0.99] = 0.99   

            # set up empty column size distribution array
            SPCdist = np.zeros((self.ntyme,)+R.shape)
            # loop throuhg time steps
            for t in range(self.ntyme):        
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
                    SPCdist[t,:] = SPCdist[t,:] + dvdr

            self.__dict__[spc+'dist'].append(SPCdist)


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
            chmin = np.argmax(dch[dch<0])
            chmin = MODIS_channels[dch<0][chmin]
            chmax = np.argmin(dch[dch>0])
            chmax = MODIS_channels[dch>0][chmax]

            SDS = 'Riso'+chmin,'Rgeo'+chmin,'Rvol'+chmin,'Riso'+chmax,'Rgeo'+chmax,'Rvol'+chmax

        if self.verbose:
            print 'opening BRDF file ',self.brdfFile
        nc = Dataset(self.brdfFile)

        for sds in SDS:
            self.__dict__[sds] = np.array(nc.variables[sds][:])

        missing_value = nc.variables[sds].missing_value
        nc.close()
        
        # Interpolate if necessary
        if chs not in MODIS_channels:
            X = np.array([chmin,chmax]).astype('int')
            for R in ['Riso','Rgeo','Rvol']:
                sds = R + chs
                self.__dict__[sds] = np.empty([self.nstations,self.ntyme])
                for s in range(self.nstations):
                    for i in range(self.ntyme):
                        Y = np.array([self.__dict__[R+chmin][s,i],self.__dict__[R+chmax][s,i]])
                        if missing_value in Y:
                            self.__dict__[sds][s,i] = missing_value
                        else:
                            f = interpolate.interp1d(X, Y)
                            self.__dict__[sds][s,i] = f([int(chs)])

            SDS = 'Riso'+chs,'Rgeo'+chs,'Rvol'+chs
        
        # Check for missing kernel weights
        Riso = self.__dict__['Riso'+chs]
        Rgeo = self.__dict__['Rgeo'+chs]
        Rvol = self.__dict__['Rvol'+chs]
        for t in range(self.ntyme):
            iGood = (Riso[:,t] != missing_value) & (Rgeo[:,t] != missing_value) & (Rvol[:,t] != missing_value)
            self.iGood[t] = self.iGood[t] & iGood
            self.nobs[t] = np.sum(self.iGood[t])

        for sds in SDS:
            self.__dict__[sds].shape = (1,1,self.nstations,self.ntyme)

        # [nkernel,nch,nobs]
        self.kernel_wt = []
        for t in range(self.ntyme):
            kernel_wt = np.append(self.__dict__['Riso'+chs][:,:,:,t],self.__dict__['Rgeo'+chs][:,:,:,t],axis=0)
            kernel_wt = np.append(kernel_wt,self.__dict__['Rvol'+chs][:,:,:,t],axis=0)
            self.kernel_wt.append(kernel_wt)

        param1 = np.array([2]*self.nstations)
        param2 = np.array([1]*self.nstations)

        param1.shape = (1,1,self.nstations)
        param2.shape = (1,1,self.nstations)

        # [nparam,nch,nobs]
        self.RTLSparam = [np.append(param1,param2,axis=0)]*self.ntyme

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
            print 'opening BRDF file ',self.brdfFile
        nc = Dataset(self.brdfFile)

        for sds in mSDS:
            self.__dict__[sds] = np.array(nc.variables[sds][:])

        missing_value = nc.variables[sds].missing_value
        nc.close()

        if self.verbose:
            print 'opening LER albedo file ',self.lerFile
        nc = Dataset(self.lerFile)

        self.__dict__[lSDS] = np.array(nc.variables[lSDS][:,:,0])
        missing_value_l = nc.variables[lSDS].missing_value
        nc.close()

        if missing_value != missing_value_l:
            I = self.__dict__[lSDS] == missing_value_l
            self.__dict__[lSDS][I] = missing_value


        # Interpolate 
        X = np.array([388,470])
        #Riso
        sds = 'Riso' + chs
        R   = 'Riso' + MODISchs
        self.__dict__[sds] = np.empty([self.nstations,self.ntyme])
        for s in range(self.nstations):
            for i in range(self.ntyme):
                Y = np.array([self.__dict__[lSDS][s,i],self.__dict__[R][s,i]])                
                if missing_value in Y:
                    self.__dict__[sds][s,i] = missing_value
                else:
                    f = interpolate.interp1d(X, Y)
                    self.__dict__[sds][s,i] = f([int(chs)])

        #Rgeo and Rvol
        for R in 'Rgeo','Rvol':
            sds = R + chs
            self.__dict__[sds] = np.empty([self.nstations,self.ntyme])
            for s in range(self.nstations):
                for i in range(self.ntyme):
                    Y = np.array([0.0,self.__dict__[R+MODISchs][s,i]])
                    if missing_value in Y:
                        self.__dict__[sds][s,i] = missing_value
                    else:
                        f = interpolate.interp1d(X, Y)
                        self.__dict__[sds][s,i] = f([int(chs)])

        SDS = 'Riso'+chs,'Rgeo'+chs,'Rvol'+chs
        
        # Check for missing kernel weights
        Riso = self.__dict__['Riso'+chs]
        Rgeo = self.__dict__['Rgeo'+chs]
        Rvol = self.__dict__['Rvol'+chs]
        for t in range(self.ntyme):
            iGood = (Riso[:,t] != missing_value) & (Rgeo[:,t] != missing_value) & (Rvol[:,t] != missing_value)
            self.iGood[t] = self.iGood[t] & iGood
            self.nobs[t] = np.sum(self.iGood[t])

        for sds in SDS:
            self.__dict__[sds].shape = (1,1,self.nstations,self.ntyme)

        # [nkernel,nch,nobs]
        self.kernel_wt = []
        for t in range(self.ntyme):
            kernel_wt = np.append(self.__dict__['Riso'+chs][:,:,:,t],self.__dict__['Rgeo'+chs][:,:,:,t],axis=0)
            kernel_wt = np.append(kernel_wt,self.__dict__['Rvol'+chs][:,:,:,t],axis=0)
            self.kernel_wt.append(kernel_wt)

        param1 = np.array([2]*self.nstations)
        param2 = np.array([1]*self.nstations)

        param1.shape = (1,1,self.nstations)
        param2.shape = (1,1,self.nstations)

        # [nparam,nch,nobs]
        self.RTLSparam = [np.append(param1,param2,axis=0)]*self.ntyme

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
            self.__dict__[sds] = np.array(nc.variables[sds][:,:,0])

        missing_value = nc.variables[sds].missing_value
        nc.close()

        sds = 'SRFLER'+chs
        # interpolate if needed
        if chs not in LER_channels:
            dch   = self.channel - LER_channels.astype('int')
            try:
                chmin = np.argmax(dch[dch<0])
                chmin = LER_channels[dch<0][chmin]
            except:
                chmin = str(LER_channels.astype('int').min())
            try:
                chmax = np.argmin(dch[dch>0])
                chmax = LER_channels[dch>0][chmax]            
            except:                
                chmax = str(LER_channels.astype('int').max())

            X = np.array([chmin,chmax]).astype('int')

            self.__dict__[sds] = np.empty([self.nstations,self.ntyme])
            for s in range(self.nstations):
                for i in range(self.ntyme):
                    Y = np.array([self.__dict__['SRFLER'+chmin][s,i],self.__dict__['SRFLER'+chmax][s,i]])
                    if missing_value in Y:                        
                        self.__dict__[sds][s,i] = missing_value
                    else:
                        try:
                            f = interpolate.interp1d(X, Y,fill_value='extrapolate')
                            self.__dict__[sds][s,i] = f([int(chs)]) 
                        except:
                            f = extrap1d(interpolate.interp1d(X, Y))
                            self.__dict__[sds][s,i] = f([int(chs)]) 


        # Check for missing values
        for t in range(self.ntyme):
            iGood = (self.__dict__[sds][:,t] != missing_value) 
            self.iGood[t] = self.iGood[t] & iGood
            self.nobs[t] = np.sum(self.iGood[t])

        # (nobs,nch)
        self.__dict__[sds].shape = (self.nstations,self.ntyme,1) 

        self.albedo = []
        for t in range(self.ntyme):
            self.albedo.append(self.__dict__[sds][:,t,:])

    #---
    def computeAE(self):
        """
        Computes angstrom exponent between 440 and 870 nm
        """
        self.ae = []

        # get 440 tau
        self.channel = 440
        self.rcFile = 'Aod_EOS.440.rc'
        pool = multiprocessing.Pool(int(multiprocessing.cpu_count()*0.5))     
        args = zip([self]*self.ntyme,range(self.ntyme))
        result = pool.map(unwrap_self_doMie,args)

        tau440 = []        
        for r in result:
            tau, ssa, g, pmom,refi,refr,vol = r
            tau = np.squeeze(tau)
            tau440.append(tau.sum(axis=0))  #(nobs)

        pool.close()
        pool.join()


        # # get 870 tau
        self.channel = 870
        self.rcFile = 'Aod_EOS.870.rc'
        pool = multiprocessing.Pool(int(multiprocessing.cpu_count()*0.5))     
        args = zip([self]*self.ntyme,range(self.ntyme))
        result = pool.map(unwrap_self_doMie,args)

        tau870 = []        
        for r in result:
            tau, ssa, g, pmom,refi,refr,vol = r
            tau = np.squeeze(tau)
            tau870.append(tau.sum(axis=0))  #(nobs)

        ratio = np.log(870./440.)
        for t870,t440 in zip(tau870,tau440):
            ae = -1.*np.log(t870/t440)/ratio
            self.ae.append(ae)

        pool.close()
        pool.join()


    #---
    def computeMie(self):
        """
        Computes aerosol optical quantities 
        """
        self.tau = []
        self.ssa = []
        self.g   = []
        self.pmom = []
        self.refi = []
        self.refr = []
        self.vol  = []

        pool = multiprocessing.Pool(int(multiprocessing.cpu_count()*0.5))     
        args = zip([self]*self.ntyme,range(self.ntyme))
        result = pool.map(unwrap_self_doMie,args)
        
        for r in result:
            tau, ssa, g, pmom,refi,refr,vol = r
            self.tau.append(tau)  #(km,nch,nobs)
            self.ssa.append(ssa)  #(km,nch,nobs)
            self.g.append(g)    #(km,nch,nobs)
            self.pmom.append(pmom)  #(km,nch,nobs,nMom,nPol)
            self.refi.append(refi) #(km,nch,nobs)
            self.refr.append(refr) #(km,nch,nobs)
            self.vol.append(vol) #(km,nch,nobs)

        pool.close()
        pool.join()

        

    def doMie(self,t):
        NAMES = self.AERNAMES + ['PS','DELP','RH','AIRDENS']
        inVars = MieVARS()
        for sds in NAMES:
            var = self.__dict__[sds]
            if len(var.shape)==3:
                inVars.__dict__[sds] = var[:,t,:]
            elif len(var.shape) ==2:
                inVars.__dict__[sds] = var[:,t]

        tau,ssa,g,pmom = getAOPvector(inVars,self.channel,
                                 vnames=self.AERNAMES,
                                 Verbose=False,
                                 rcfile=self.rcFile,
                                 nMom=self.nMom)

        vol, area, refr, refi, reff = getAOPint(inVars,self.channel,
                                 vnames=self.AERNAMES,
                                 Verbose=False,
                                 rcfile=self.rcFile)

        return tau,ssa,g,pmom,refr,refi,vol

    # --
    def runExt(self):
        """
        run ext_sampler.py 
        """
        if self.verbose: 
            print 'running ext_sampler on file',self.inFile.replace('%col',self.extcol)

        Options =     " --input=" + self.inFile.replace('%col',self.extcol)      + \
                      " --output=" + self.outFileext.replace('%col',self.extcol)       + \
                      " --rc=" + self.rcFile      + \
                      " --format=NETCDF4_CLASSIC"      + \
                      " --channel=%d" %self.channel + \
                      " --intensive"   +\
                      " --stn" +\
                      " --oc" +\
                      " --bc" +\
                      " --ss" +\
                      " --su" +\
                      " --du"    
                      

        if not os.path.exists(os.path.dirname(self.outFileext)):
            os.makedirs(os.path.dirname(self.outFileext))

        cmd = 'ext_stn_vlidort_sampler.py {} '.format(Options)  
        print cmd
        if os.system(cmd):
            raise ValueError, "ext_stn_vlidort_sampler.py failed for %s "%(self.inFile.replace('%col',self.extcol))       


    def runALMUCANTAR(self):
        """
        Calls VLIDORT for almuncantar scan
        """

        # index of stations with good obs
        tymes     = np.arange(self.ntyme)[self.nobs>0]

        self.I = []
        self.Q = []
        self.U = []
        self.reflectance = []
        self.surf_reflectance = []
        self.ROT = []
        for t in tymes: 
            tau = self.tau[t][:,:,self.iGood[t]]
            ssa = self.ssa[t][:,:,self.iGood[t]]
            pmom = self.pmom[t][:,:,self.iGood[t],:,:]
            pe   = self.pe[t][:,self.iGood[t]]
            ze   = self.ze[t][:,self.iGood[t]]
            te   = self.te[t][:,self.iGood[t]]

            # Initiate output arrays
            nlev    = tau.shape[0]
            self.I_ = np.ones([self.nstations,self.nal])*MISSING
            self.Q_ = np.ones([self.nstations,self.nal])*MISSING
            self.U_ = np.ones([self.nstations,self.nal])*MISSING
            self.reflectance_ = np.ones([self.nstations,self.nal])*MISSING
            self.surf_reflectance_ = np.ones([self.nstations,self.nal])*MISSING
            self.ROT_ = np.ones([self.nstations,nlev])*MISSING
        
            # Get VLIDORT wrapper name from dictionary
            vlidortWrapper = WrapperFuncs[self.albedoType]


            # do almuncantar
            for ai, a in enumerate(self.al_angles):         
                # Solar Geometry
                sza  = self.SZA[t][self.iGood[t]]
                saa  = self.SAA[t][self.iGood[t]]

                #viewing geometry
                vza  = sza                
                vaa  = saa + a
                iGood = vaa < 0
                vaa[iGood] = 360.0 + vaa[iGood]
                iGood = vaa >= 360.0
                vaa[iGood] = vaa[iGood] - 360.00

                # define RAA according to photon travel direction
                saa = saa + 180.0
                iGood = saa >= 360.0
                saa[iGood] = saa[iGood] - 360.0

                raa = np.abs(vaa - saa)

                # run VLIDORT                
                # get args list for each surface model
                if self.albedoType == 'MODIS_BRDF':
                    kernel_wt = self.kernel_wt[t][:,:,self.iGood[t]]
                    param     = self.RTLSparam[t][:,:,self.iGood[t]]                
                    
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
                    kernel_wt = self.kernel_wt[t][:,:,self.iGood[t]]
                    RTLSparam = self.RTLSparam[t][:,:,self.iGood[t]] 
                    RTLSparam = np.append(RTLSparam,np.zeros([1,1,self.nobs[t]]),axis=0) 

                    # For BPDF
                    BPDFparam = self.BPDFparam[t][:,:,self.iGood[t]]

                    # Loop through one by one
                    # Some land covers do not have polarization (i.e. urban)
                    I = np.zeros([self.nobs[t],1])
                    Q = np.zeros([self.nobs[t],1])
                    U = np.zeros([self.nobs[t],1])
                    reflectance = np.zeros([self.nobs[t],1])
                    surf_reflectance = np.zeros([self.nobs[t],1])
                    ROT = np.zeros([nlev,self.nobs[t],1])

                    for p in range(self.nobs[t]):
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
                    albedo = self.albedo[t][self.iGood[t],:]
                    
                    args = [self.channel,tau, ssa, pmom, 
                            pe, ze, te, 
                            albedo, 
                            sza, raa, vza, 
                            MISSING,
                            self.verbose]

                    # Call VLIDORT wrapper function
                    I, reflectance, ROT, Q, U, rc = vlidortWrapper(*args)  
                    surf_reflectance = albedo

                # Store values in initialized arrays
                if ai == 0:
                    self.ROT_[self.iGood[t],:] = np.squeeze(ROT).T

                self.I_[self.iGood[t],ai]                = np.squeeze(I)
                self.reflectance_[self.iGood[t],ai]      = np.squeeze(reflectance)
                self.surf_reflectance_[self.iGood[t],ai] = np.squeeze(surf_reflectance)
                self.Q_[self.iGood[t],ai]                = np.squeeze(Q)
                self.U_[self.iGood[t],ai]                = np.squeeze(U) 

            #Store full almucantar for writing later
            self.I.append(self.I_)
            self.Q.append(self.Q_)
            self.U.append(self.U_)
            self.reflectance.append(self.reflectance_)
            self.surf_reflectance.append(self.surf_reflectance_)
            self.ROT.append(self.ROT_)

        self.writeNCal()

    def runPRINCIPLE(self):
        """
        Calls VLIDORT for principle place scan
        """

        # index of stations with good obs
        tymes     = np.arange(self.ntyme)[self.nobs>0]

        self.I = []
        self.Q = []
        self.U = []
        self.reflectance = []
        self.surf_reflectance = []
        self.ROT = []
        for t in tymes: 
            tau = self.tau[t][:,:,self.iGood[t]]
            ssa = self.ssa[t][:,:,self.iGood[t]]
            pmom = self.pmom[t][:,:,self.iGood[t],:,:]
            pe   = self.pe[t][:,self.iGood[t]]
            ze   = self.ze[t][:,self.iGood[t]]
            te   = self.te[t][:,self.iGood[t]]

            # Initiate output arrays
            nlev    = tau.shape[0]
            self.I_ = np.ones([self.nstations,self.npp])*MISSING
            self.Q_ = np.ones([self.nstations,self.npp])*MISSING
            self.U_ = np.ones([self.nstations,self.npp])*MISSING
            self.reflectance_ = np.ones([self.nstations,self.npp])*MISSING
            self.surf_reflectance_ = np.ones([self.nstations,self.npp])*MISSING
            self.ROT_ = np.ones([self.nstations,nlev])*MISSING
        
            # Get VLIDORT wrapper name from dictionary
            vlidortWrapper = WrapperFuncs[self.albedoType]


            # do principle plane
            for pi, pp in enumerate(self.pp_angles):         
                # Solar Geometry
                sza  = self.SZA[t][self.iGood[t]]
                saa  = self.SAA[t][self.iGood[t]]

                #viewing geometry
                vza  = sza - pp  

                raa  = np.zeros(saa.shape)

                # if pointing towards the sun
                iGood = vza >= 0
                raa[iGood] = 180.0

                #if poinitng away from the sun
                iGood = vza < 0
                raa[iGood] = 0.0

                # make all vza's positive
                vza[iGood] = np.abs(vza[iGood])

                # Limit viewing angles to less than 80
                # also limits to pointing above the horizon
                iGood = vza < 80
                nobs = np.sum(iGood)

                if nobs >0:
                    # run VLIDORT                
                    # get args list for each surface model
                    if self.albedoType == 'MODIS_BRDF':
                        kernel_wt = self.kernel_wt[t][:,:,self.iGood[t]]
                        param     = self.RTLSparam[t][:,:,self.iGood[t]]                
                        
                        args = [self.channel,tau[:,:,iGood], ssa[:,:,iGood], pmom[:,:,iGood,:,:], 
                                pe[:,iGood], ze[:,iGood], te[:,iGood], 
                                kernel_wt, param, 
                                sza[iGood], raa[iGood], vza[iGood], 
                                MISSING,
                                self.verbose]

                        # Call VLIDORT wrapper function
                        I, reflectance, ROT, surf_reflectance, Q, U, rc = vlidortWrapper(*args)                        
                        
                    elif self.albedoType == 'MODIS_BRDF_BPDF':
                        # For albedo
                        kernel_wt = self.kernel_wt[t][:,:,self.iGood[t]]
                        RTLSparam = self.RTLSparam[t][:,:,self.iGood[t]] 
                        RTLSparam = np.append(RTLSparam,np.zeros([1,1,self.nobs[t]]),axis=0) 

                        # For BPDF
                        BPDFparam = self.BPDFparam[t][:,:,self.iGood[t]]

                        # Loop through one by one
                        # Some land covers do not have polarization (i.e. urban)
                        I = np.zeros([nobs,1])
                        Q = np.zeros([nobs,1])
                        U = np.zeros([nobs,1])
                        reflectance = np.zeros([nobs,1])
                        surf_reflectance = np.zeros([nobs,1])
                        ROT = np.zeros([nlev,nobs,1])

                        index = np.arange(self.nobs[t])
                        index = index[iGood]

                        for pindex,p in enumerate(index):
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
                    
                            I[pindex:pindex+1,:] = I_
                            Q[pindex:pindex+1,:] = Q_
                            U[pindex:pindex+1,:] = U_
                            reflectance[pindex:pindex+1,:] = reflectance_
                            surf_reflectance[pindex:pindex+1,:] = surf_reflectance_
                            ROT[:,pindex:pindex+1,:] = ROT_
                    
                    elif self.albedoType == 'LAMBERTIAN':
                        albedo = self.albedo[t][self.iGood[t],:]
                        
                        args = [self.channel,tau[:,:,iGood], ssa[:,:,iGood], pmom[:,:,iGood,:,:], 
                                pe[:,iGood], ze[:,iGood], te[:,iGood], 
                                albedo[iGood,:], 
                                sza[iGood], raa[iGood], vza[iGood], 
                                MISSING,
                                self.verbose]

                        # Call VLIDORT wrapper function
                        I, reflectance, ROT, Q, U, rc = vlidortWrapper(*args)  
                        surf_reflectance = albedo[iGood,:]

                    # Store values in initialized arrays
                    II = np.arange(self.nstations)
                    II = II[self.iGood[t]]
                    II = II[iGood]
                    if pi == 0:
                        self.ROT_[II,:] = np.squeeze(ROT).T

                    self.I_[II,pi]                = np.squeeze(I)
                    self.reflectance_[II,pi]      = np.squeeze(reflectance)
                    self.surf_reflectance_[II,pi] = np.squeeze(surf_reflectance)
                    self.Q_[II,pi]                = np.squeeze(Q)
                    self.U_[II,pi]                = np.squeeze(U) 

            #Store full principle plane for writing later
            self.I.append(self.I_)
            self.Q.append(self.Q_)
            self.U.append(self.U_)
            self.reflectance.append(self.reflectance_)
            self.surf_reflectance.append(self.surf_reflectance_)
            self.ROT.append(self.ROT_)

        self.writeNCpp()



        
    #---
    def writeNCal (self,zlib=True):
        """
        Write a NetCDF file vlidort output
        """
        km = 72

        if not os.path.exists(os.path.dirname(self.outFileal)):
            os.makedirs(os.path.dirname(self.outFileal))

        # Open NC file
        # ------------
        nc = Dataset(self.outFileal,'w',format='NETCDF4_CLASSIC')

        # Set global attributes
        # ---------------------
        nc.title = 'VLIDORT Simulation of GEOS-5 almucantar scan'
        nc.institution = 'NASA/Goddard Space Flight Center'
        nc.source = 'Global Model and Assimilation Office'
        nc.history = 'VLIDORT simulation run on sampled GEOS-5'
        nc.references = 'n/a'
        nc.note = 'Sun normalized output.  Solar flux set equal to 1.'
        nc.contact = 'Patricia Castellanos <patricia.castellanos@nasa.gov>'
        nc.Conventions = 'CF'
        nc.inFile = self.inFile
     
        # Create dimensions
        # -----------------
        nt = nc.createDimension('time',len(self.tyme))
        x  = nc.createDimension('x',1)
        y  = nc.createDimension('y',1)
        ns = nc.createDimension('station',self.nstations)
        na = nc.createDimension('angle',self.nal)        
        nz = nc.createDimension('lev',km)
        ls = nc.createDimension('ls',19)

        # Coordinate variables
        # --------------------
        col = 'aer_Nv'
        if self.verbose: 
            print 'opening file',self.inFile.replace('%col',col)
        nctrj       = Dataset(self.inFile.replace('%col',col))     
        _copyVar(nctrj,nc,u'time', dtype='i4',zlib=False,verbose=self.verbose)   
        _copyVar(nctrj,nc,u'x',dtype='f4',zlib=False,verbose=self.verbose)
        _copyVar(nctrj,nc,u'y',dtype='f4',zlib=False,verbose=self.verbose)  
        _copyVar(nctrj,nc,u'station',dtype='f4',zlib=False,verbose=self.verbose)   

        aa = nc.createVariable('angle','f4',('angle',),zlib=False)
        aa.long_name     = 'azimuth angle relative to sun'
        aa.units         = 'degrees'
        aa[:]            = np.array(self.al_angles)

        _copyVar(nctrj,nc,u'lev',dtype='f4',zlib=False,verbose=self.verbose)       
        _copyVar(nctrj,nc,u'stnLon',dtype='f4',zlib=False,verbose=self.verbose)
        _copyVar(nctrj,nc,u'stnLat',dtype='f4',zlib=False,verbose=self.verbose)
        _copyVar(nctrj,nc,u'stnName',dtype='S1',zlib=False,verbose=self.verbose)
        _copyVar(nctrj,nc,u'isotime', dtype='S1',zlib=False,verbose=self.verbose)

        nctrj.close()

        # Write VLIDORT Outputs
        # ---------------------
        # dummy array to hold data
        b = np.empty([self.nstations,self.nal,self.ntyme])


        ref = nc.createVariable('I_normalized','f4',('station','angle','time',),zlib=zlib,fill_value=MISSING)
        ref.standard_name = '%.2f nm BOA normalized I' %self.channel
        ref.long_name     = '%.2f nm normalized intensity at the bottom of the atmosphere (pi*I/cosSZA*SOLAR_FLUX)' %self.channel
        ref.missing_value = MISSING
        ref.units         = "None"
        for t in range(self.ntyme):
            b[:,:,t] = self.reflectance[t]
        ref[:]            = b

        i = nc.createVariable('I','f4',('station','angle','time',),zlib=zlib,fill_value=MISSING)
        i.standard_name = '%.2f nm BOA I' %self.channel
        i.long_name     = '%.2f nm intensity at the bottom of the atmosphere' %self.channel
        i.missing_value = MISSING
        i.units         = "W m-2 sr-1 nm-1"
        for t in range(self.ntyme):
            b[:,:,t] = self.I[t]
        i[:]            = b

        q = nc.createVariable('Q','f4',('station','angle','time',),zlib=zlib,fill_value=MISSING)
        q.standard_name = '%.2f nm BOA Q' %self.channel
        q.long_name     = '%.2f nm Q-component of the stokes vector at the bottom of the atmopshere' %self.channel
        q.missing_value = MISSING
        q.units         = "W m-2 sr-1 nm-1"
        for t in range(self.ntyme):
            b[:,:,t] = self.Q[t]
        q[:]            = b

        u = nc.createVariable('U','f4',('station','angle','time',),zlib=zlib,fill_value=MISSING)
        u.standard_name = '%.2f nm BOA U' %self.channel
        u.long_name     = '%.2f nm U-component of the stokes vector at the bottom of the atmopshere' %self.channel
        u.missing_value = MISSING
        u.units         = "W m-2 sr-1 nm-1"
        for t in range(self.ntyme):
            b[:,:,t] = self.U[t]        
        u[:]            = b

        sref = nc.createVariable('surf_reflectance','f4',('station','angle','time',),zlib=zlib,fill_value=MISSING)
        sref.standard_name = '%.2f nm Surface Reflectance' %self.channel
        sref.long_name     = '%.2f nm Bi-Directional Surface Reflectance' %self.channel
        sref.missing_value = MISSING
        sref.units         = "None"
        for t in range(self.ntyme):
            b[:,:,t] = self.surf_reflectance[t]
        sref[:]            = b

        # SSA
        ssa = nc.createVariable('ssa','f4',('station','time',),zlib=zlib,fill_value=MISSING)
        ssa.standard_name = '%.2f nm aerosol SSA' %self.channel
        ssa.long_name     = '%.2f nm aerosol single scattering albedo' %self.channel
        ssa.missing_value = MISSING
        ssa.units         = "None"
        b = np.empty([self.nstations,self.ntyme])
        for t in range(self.ntyme):
            ssatemp  = np.squeeze(self.ssa[t])
            tau      = np.squeeze(self.tau[t])
            b[:,t] = (ssatemp*tau).sum(axis=0)/tau.sum(axis=0)
        ssa[:]            = b

        # AOD
        aod = nc.createVariable('aod','f4',('station','time',),zlib=zlib,fill_value=MISSING)
        aod.standard_name = '%.2f nm AOD' %self.channel
        aod.long_name     = '%.2f nm aerosol optical depth' %self.channel
        aod.missing_value = MISSING
        aod.units         = "None"
        b = np.empty([self.nstations,self.ntyme])
        for t in range(self.ntyme):
            tau      = np.squeeze(self.tau[t])
            b[:,t]   = tau.sum(axis=0)
        aod[:]            = b

        # Refractive Index
        # Get column refractive index
        b = np.empty([self.nstations,self.ntyme])
        for t in range(self.ntyme):    
            vol = np.squeeze(self.vol[t])
            refi = np.squeeze(self.refi[t])
            b[:,t] = (vol*refi).sum(axis=0)/vol.sum(axis=0)

        
        refi = nc.createVariable('refi','f4',('station','time',),zlib=zlib,fill_value=MISSING)
        refi.standard_name = '%.2f aerosol imaginary refractive index' %self.channel
        refi.long_name     = '%.2f aerosol imaginary refractive index' %self.channel
        refi.missing_value = MISSING
        refi.units         = "None"
        refi[:]            = b

        b = np.empty([self.nstations,self.ntyme])
        for t in range(self.ntyme):    
            vol = np.squeeze(self.vol[t])
            refr = np.squeeze(self.refr[t])
            b[:,t] = (vol*refr).sum(axis=0)/vol.sum(axis=0)

        
        refr = nc.createVariable('refr','f4',('station','time',),zlib=zlib,fill_value=MISSING)
        refr.standard_name = '%.2f aerosol real refractive index' %self.channel
        refr.long_name     = '%.2f aerosol real refractive index' %self.channel
        refr.missing_value = MISSING
        refr.units         = "None"
        refr[:]            = b




        if "MODIS_BRDF" in self.albedoType:
            # BRDF Params
            chs = str(int(self.channel))

            riso = nc.createVariable('riso','f4',('station','time',),zlib=zlib,fill_value=MISSING)
            riso.standard_name = '%.2f nm RISO' %self.channel
            riso.long_name     = '%.2f nm isotropic BRDF kernel' %self.channel
            riso.missing_value = MISSING
            riso.units         = "None"
            riso[:]            = self.__dict__['Riso'+chs][0,0,:,:]


            rgeo = nc.createVariable('rgeo','f4',('station','time',),zlib=zlib,fill_value=MISSING)
            rgeo.standard_name = '%.2f nm RGEO' %self.channel
            rgeo.long_name     = '%.2f nm geometric BRDF kernel' %self.channel
            rgeo.missing_value = MISSING
            rgeo.units         = "None"
            rgeo[:]            = self.__dict__['Rgeo'+chs][0,0,:,:]

            rvol = nc.createVariable('rvol','f4',('station','time',),zlib=zlib,fill_value=MISSING)
            rvol.standard_name = '%.2f nm RVOL' %self.channel
            rvol.long_name     = '%.2f nm volumetric BRDF kernel' %self.channel
            rvol.missing_value = MISSING
            rvol.units         = "None"
            rvol[:]            = self.__dict__['Rvol'+chs][0,0,:,:]




        # Vertical Profile Variables
        # dummy array to hold data
        b = np.empty([self.nstations,self.ntyme,km])

        rot = nc.createVariable('ROT','f4',('station','time','lev',),zlib=zlib,fill_value=MISSING)
        rot.long_name = '%.2f nm Rayleigh Optical Thickness' %self.channel
        rot.missing_value = MISSING
        rot.units         = "None"
        for t in range(self.ntyme):
            b[:,t,:] = self.ROT[t]
        rot[:]            = b


        # Close the file
        # --------------
        nc.close()

        if self.verbose:
            print " <> wrote %s"%(self.outFileal)
    

    #---
    def writeNCpp (self,zlib=True):
        """
        Write a NetCDF file vlidort output
        """
        km = 72

        if not os.path.exists(os.path.dirname(self.outFilepp)):
            os.makedirs(os.path.dirname(self.outFilepp))

        # Open NC file
        # ------------
        nc = Dataset(self.outFilepp,'w',format='NETCDF4_CLASSIC')

        # Set global attributes
        # ---------------------
        nc.title = 'VLIDORT Simulation of GEOS-5 principle plane scan'
        nc.institution = 'NASA/Goddard Space Flight Center'
        nc.source = 'Global Model and Assimilation Office'
        nc.history = 'VLIDORT simulation run on sampled GEOS-5'
        nc.references = 'n/a'
        nc.note = 'Sun normalized output.  Solar flux set equal to 1.'
        nc.contact = 'Patricia Castellanos <patricia.castellanos@nasa.gov>'
        nc.Conventions = 'CF'
        nc.inFile = self.inFile
     
        # Create dimensions
        # -----------------
        nt = nc.createDimension('time',len(self.tyme))
        x  = nc.createDimension('x',1)
        y  = nc.createDimension('y',1)
        ns = nc.createDimension('station',self.nstations)
        na = nc.createDimension('angle',self.npp)        
        nz = nc.createDimension('lev',km)
        ls = nc.createDimension('ls',19)

        # Coordinate variables
        # --------------------
        col = 'aer_Nv'
        if self.verbose: 
            print 'opening file',self.inFile.replace('%col',col)
        nctrj       = Dataset(self.inFile.replace('%col',col))     
        _copyVar(nctrj,nc,u'time', dtype='i4',zlib=False,verbose=self.verbose)   
        _copyVar(nctrj,nc,u'x',dtype='f4',zlib=False,verbose=self.verbose)
        _copyVar(nctrj,nc,u'y',dtype='f4',zlib=False,verbose=self.verbose)  
        _copyVar(nctrj,nc,u'station',dtype='f4',zlib=False,verbose=self.verbose)   

        aa = nc.createVariable('angle','f4',('angle',),zlib=False)
        aa.long_name     = 'Scattering Angle from sun (negative is below the sun) '
        aa.units         = 'degrees'
        aa[:]            = np.array(self.pp_angles)

        _copyVar(nctrj,nc,u'lev',dtype='f4',zlib=False,verbose=self.verbose)       
        _copyVar(nctrj,nc,u'stnLon',dtype='f4',zlib=False,verbose=self.verbose)
        _copyVar(nctrj,nc,u'stnLat',dtype='f4',zlib=False,verbose=self.verbose)
        _copyVar(nctrj,nc,u'stnName',dtype='S1',zlib=False,verbose=self.verbose)
        _copyVar(nctrj,nc,u'isotime', dtype='S1',zlib=False,verbose=self.verbose)

        nctrj.close()

        # Write VLIDORT Outputs
        # ---------------------
        # dummy array to hold data
        b = np.empty([self.nstations,self.npp,self.ntyme])


        ref = nc.createVariable('I_normalized','f4',('station','angle','time',),zlib=zlib,fill_value=MISSING)
        ref.standard_name = '%.2f nm BOA normalized I' %self.channel
        ref.long_name     = '%.2f nm normalized intensity at the bottom of the atmosphere (pi*I/cosSZA*SOLAR_FLUX)' %self.channel
        ref.missing_value = MISSING
        ref.units         = "None"
        for t in range(self.ntyme):
            b[:,:,t] = self.reflectance[t]
        ref[:]            = b

        i = nc.createVariable('I','f4',('station','angle','time',),zlib=zlib,fill_value=MISSING)
        i.standard_name = '%.2f nm BOA I' %self.channel
        i.long_name     = '%.2f nm intensity at the bottom of the atmosphere' %self.channel
        i.missing_value = MISSING
        i.units         = "W m-2 sr-1 nm-1"
        for t in range(self.ntyme):
            b[:,:,t] = self.I[t]
        i[:]            = b

        q = nc.createVariable('Q','f4',('station','angle','time',),zlib=zlib,fill_value=MISSING)
        q.standard_name = '%.2f nm BOA Q' %self.channel
        q.long_name     = '%.2f nm Q-component of the stokes vector at the bottom of the atmopshere' %self.channel
        q.missing_value = MISSING
        q.units         = "W m-2 sr-1 nm-1"
        for t in range(self.ntyme):
            b[:,:,t] = self.Q[t]
        q[:]            = b

        u = nc.createVariable('U','f4',('station','angle','time',),zlib=zlib,fill_value=MISSING)
        u.standard_name = '%.2f nm BOA U' %self.channel
        u.long_name     = '%.2f nm U-component of the stokes vector at the bottom of the atmopshere' %self.channel
        u.missing_value = MISSING
        u.units         = "W m-2 sr-1 nm-1"
        for t in range(self.ntyme):
            b[:,:,t] = self.U[t]        
        u[:]            = b

        sref = nc.createVariable('surf_reflectance','f4',('station','angle','time',),zlib=zlib,fill_value=MISSING)
        sref.standard_name = '%.2f nm Surface Reflectance' %self.channel
        sref.long_name     = '%.2f nm Bi-Directional Surface Reflectance' %self.channel
        sref.missing_value = MISSING
        sref.units         = "None"
        for t in range(self.ntyme):
            b[:,:,t] = self.surf_reflectance[t]
        sref[:]            = b

        ssa = nc.createVariable('ssa','f4',('station','time',),zlib=zlib,fill_value=MISSING)
        ssa.standard_name = '%.2f nm aerosol single scattering albedo' %self.channel
        ssa.long_name     = '%.2f nm aerosol single scattering albedo' %self.channel
        ssa.missing_value = MISSING
        ssa.units         = "None"
        b = np.empty([self.nstations,self.ntyme])
        for t in range(self.ntyme):
            ssatemp  = np.squeeze(self.ssa[t])
            tau      = np.squeeze(self.tau[t])
            b[:,t] = (ssatemp*tau).sum(axis=0)/tau.sum(axis=0)
        ssa[:]            = b

        aod = nc.createVariable('aod','f4',('station','time',),zlib=zlib,fill_value=MISSING)
        aod.standard_name = '%.2f nm AOD' %self.channel
        aod.long_name     = '%.2f nm aerosol optical depth' %self.channel
        aod.missing_value = MISSING
        aod.units         = "None"
        b = np.empty([self.nstations,self.ntyme])
        for t in range(self.ntyme):
            tau      = np.squeeze(self.tau[t])
            b[:,t]   = tau.sum(axis=0)
        aod[:]            = b

        # Refractive Index
        # Get column refractive index
        b = np.empty([self.nstations,self.ntyme])
        for t in range(self.ntyme):    
            vol = np.squeeze(self.vol[t])
            refi = np.squeeze(self.refi[t])
            b[:,t] = (vol*refi).sum(axis=0)/vol.sum(axis=0)

        
        refi = nc.createVariable('refi','f4',('station','time',),zlib=zlib,fill_value=MISSING)
        refi.standard_name = '%.2f aerosol imaginary refractive index' %self.channel
        refi.long_name     = '%.2f aerosol imaginary refractive index' %self.channel
        refi.missing_value = MISSING
        refi.units         = "None"
        refi[:]            = b

        b = np.empty([self.nstations,self.ntyme])
        for t in range(self.ntyme):    
            vol = np.squeeze(self.vol[t])
            refr = np.squeeze(self.refi[t])
            b[:,t] = (vol*refr).sum(axis=0)/vol.sum(axis=0)

        
        refr = nc.createVariable('refr','f4',('station','time',),zlib=zlib,fill_value=MISSING)
        refr.standard_name = '%.2f aerosol real refractive index' %self.channel
        refr.long_name     = '%.2f aerosol real refractive index' %self.channel
        refr.missing_value = MISSING
        refr.units         = "None"
        refr[:]            = b


        if "MODIS_BRDF" in self.albedoType:
            # BRDF Params
            chs = str(int(self.channel))

            riso = nc.createVariable('riso','f4',('station','time',),zlib=zlib,fill_value=MISSING)
            riso.standard_name = '%.2f nm RISO' %self.channel
            riso.long_name     = '%.2f nm isotropic BRDF kernel' %self.channel
            riso.missing_value = MISSING
            riso.units         = "None"
            riso[:]            = self.__dict__['Riso'+chs][0,0,:,:]


            rgeo = nc.createVariable('rgeo','f4',('station','time',),zlib=zlib,fill_value=MISSING)
            rgeo.standard_name = '%.2f nm RGEO' %self.channel
            rgeo.long_name     = '%.2f nm geometric BRDF kernel' %self.channel
            rgeo.missing_value = MISSING
            rgeo.units         = "None"
            rgeo[:]            = self.__dict__['Rgeo'+chs][0,0,:,:]

            rvol = nc.createVariable('rvol','f4',('station','time',),zlib=zlib,fill_value=MISSING)
            rvol.standard_name = '%.2f nm RVOL' %self.channel
            rvol.long_name     = '%.2f nm volumetric BRDF kernel' %self.channel
            rvol.missing_value = MISSING
            rvol.units         = "None"
            rvol[:]            = self.__dict__['Rvol'+chs][0,0,:,:]


        # Vertical Profile Variables
        # dummy array to hold data
        b = np.empty([self.nstations,self.ntyme,km])

        rot = nc.createVariable('ROT','f4',('station','time','lev',),zlib=zlib,fill_value=MISSING)
        rot.long_name = '%.2f nm Rayleigh Optical Thickness' %self.channel
        rot.missing_value = MISSING
        rot.units         = "None"
        for t in range(self.ntyme):
            b[:,t,:] = self.ROT[t]
        rot[:]            = b

        # Close the file
        # --------------
        nc.close()

        if self.verbose:
            print " <> wrote %s"%(self.outFilepp)


    #---
    def writeNCadd (self,zlib=True):
        """
        Write a NetCDF file of additional data 
        this is all wavelength independent stuff
        """
        km = 72

        if not os.path.exists(os.path.dirname(self.outFileadd)):
            os.makedirs(os.path.dirname(self.outFileadd))

        # Open NC file
        # ------------
        nc = Dataset(self.outFileadd,'w',format='NETCDF4_CLASSIC')

        # Set global attributes
        # ---------------------
        nc.title = 'VLIDORT Simulation of GEOS-5 - Additinal Data'
        nc.institution = 'NASA/Goddard Space Flight Center'
        nc.source = 'Global Model and Assimilation Office'
        nc.history = 'VLIDORT simulation run on sampled GEOS-5'
        nc.references = 'n/a'
        nc.contact = 'Patricia Castellanos <patricia.castellanos@nasa.gov>'
        nc.Conventions = 'CF'
        nc.inFile = self.inFile
     
        # Create dimensions
        # -----------------
        nt = nc.createDimension('time',len(self.tyme))
        x  = nc.createDimension('x',1)
        y  = nc.createDimension('y',1)
        ns = nc.createDimension('station',self.nstations)
        nz = nc.createDimension('lev',km)
        ne = nc.createDimension('leve',km+1)
        ls = nc.createDimension('ls',19)
        nr = nc.createDimension('radius',len(self.aeronet_r))

        # Coordinate variables
        # --------------------
        col = 'aer_Nv'
        if self.verbose: 
            print 'opening file',self.inFile.replace('%col',col)
        nctrj       = Dataset(self.inFile.replace('%col',col))     
        _copyVar(nctrj,nc,u'time', dtype='i4',zlib=False,verbose=self.verbose)   
        _copyVar(nctrj,nc,u'x',dtype='f4',zlib=False,verbose=self.verbose)
        _copyVar(nctrj,nc,u'y',dtype='f4',zlib=False,verbose=self.verbose)  
        _copyVar(nctrj,nc,u'station',dtype='f4',zlib=False,verbose=self.verbose)   
        _copyVar(nctrj,nc,u'lev',dtype='f4',zlib=False,verbose=self.verbose)       

        leve = nc.createVariable('leve','f4',('leve',),zlib=False)
        leve.long_name   = 'Vertical Level Edge'
        leve.units        = 'layer'
        leve[:]          = np.arange(km+1)

        rad = nc.createVariable('radius','f4',('radius',),zlib=False)
        rad.long_name   = 'aerosol radius'
        rad.units       = 'microns'
        rad[:]          = self.aeronet_r*1e6


        _copyVar(nctrj,nc,u'stnLon',dtype='f4',zlib=False,verbose=self.verbose)
        _copyVar(nctrj,nc,u'stnLat',dtype='f4',zlib=False,verbose=self.verbose)
        _copyVar(nctrj,nc,u'stnName',dtype='S1',zlib=False,verbose=self.verbose)
        _copyVar(nctrj,nc,u'isotime', dtype='S1',zlib=False,verbose=self.verbose)

        nctrj.close()

        # Sun angles
        sza = nc.createVariable('solar_zenith','f4',('station','time',),zlib=zlib)
        sza.long_name     = "solar zenith angle (SZA)"
        sza.missing_value = MISSING
        sza.units         = "degrees"
        sza[:]            = np.array(self.SZA).T

        saa = nc.createVariable('solar_azimuth','f4',('station','time',),zlib=zlib)
        saa.long_name     = "solar azimuth angle (SAA)"
        saa.missing_value = MISSING
        saa.units         = "degrees clockwise from North"
        saa[:]            = np.array(self.SAA).T

        ae = nc.createVariable('angstrom_exponent','f4',('station','time',),zlib=zlib)
        ae.long_name     = "aerosol angstrom exponent between 440 and 870 nm"
        ae.missing_value = MISSING
        ae.units         = "None"
        ae[:]            = np.array(self.ae).T


        # Vertical Profile Variables
        # dummy array to hold data
        b = np.empty([self.nstations,self.ntyme,km+1])

        te = nc.createVariable('temperature','f4',('station','time','leve',),zlib=zlib,fill_value=MISSING)
        te.long_name = 'Temperature at Layer Edge' 
        te.missing_value = MISSING
        te.units         = "K"
        for t in range(self.ntyme):
            b[:,t,:] = self.te[t].T
        te[:]            = b

        pe = nc.createVariable('pressure','f4',('station','time','leve',),zlib=zlib,fill_value=MISSING)
        pe.long_name = 'Layer Edge Pressure' 
        pe.missing_value = MISSING
        pe.units         = "Pa"
        for t in range(self.ntyme):
            b[:,t,:] = self.pe[t].T
        pe[:]            = b

        ze = nc.createVariable('altitude','f4',('station','time','leve',),zlib=zlib,fill_value=MISSING)
        ze.long_name = 'Layer Edge Height Above Surface' 
        ze.missing_value = MISSING
        ze.units         = "m"
        for t in range(self.ntyme):
            b[:,t,:] = self.ze[t].T
        ze[:]            = b

        # Aerosol size distribution
        dist = nc.createVariable('aero_dist','f4',('station','time','radius',),zlib=zlib,fill_value=MISSING)
        dist.long_name     = 'aerosol size distribution (dV/dlnr)' 
        dist.missing_value = MISSING
        dist.units         = "microns^3/microns^2"
        dist[:]            = np.array(self.AERdist)


        # Close the file
        # --------------
        nc.close()

        if self.verbose:
            print " <> wrote %s"%(self.outFileadd)



    
#------------------------------------ M A I N ------------------------------------

if __name__ == "__main__":
    date     = datetime(2006,01,01,00)
    nymd     = str(date.date()).replace('-','')
    hour     = str(date.hour).zfill(2)
    format   = 'NETCDF4_CLASSIC'

    inDir        = '/nobackup/3/pcastell/STN_VLIDORT/AERONET/LevelB/Y{}/M{}'.format(date.year,str(date.month).zfill(2))
    inFile       = '{}/aeronet-g5nr.lb2.%col.{}_{}z.nc4'.format(inDir,nymd,hour)
    invDir       = '/nobackup/3/pcastell/STN_VLIDORT/AERONET/LevelB/invariant'
    invFile      = '{}/aeronet-g5nr.lb2.%col.nc4'.format(invDir)
    brdfDir      = '/nobackup/3/pcastell/STN_VLIDORT/AERONET/LevelB/surface/BRDF/MCD43C1/006/Y{}/M{}'.format(date.year,str(date.month).zfill(2))
    brdfFile     = '{}/aeronet-g5nr.lb2.brdf.{}_{}z.nc4'.format(brdfDir,nymd,hour)
    ndviDir      = '/nobackup/3/pcastell/STN_VLIDORT/AERONET/LevelB/surface/BPDF/NDVI/MYD13C2/006/Y{}/M{}'.format(date.year,str(date.month).zfill(2))
    ndviFile     = '{}/aeronet-g5nr.lb2.ndvi.{}_{}z.nc4'.format(ndviDir,nymd,hour)
    lcDir        = '/nobackup/3/pcastell/STN_VLIDORT/AERONET/LevelB/surface/BPDF/LAND_COVER/MCD12C1/051/Y{}/M{}'.format(date.year,str(date.month).zfill(2))
    lcFile       = '{}/aeronet-g5nr.lb2.land_cover.{}_{}z.nc4'.format(lcDir,nymd,hour)    
    albedoType   = 'MODIS_BRDF_BPDF'

    channel  = 470
    chd      = get_chd(channel)
    outDir    = '/nobackup/3/pcastell/STN_VLIDORT/AERONET/LevelC/Y{}/M{}'.format(date.year,str(date.month).zfill(2))
    outFileal   = '{}/aeronet-g5nr.al.vlidort.vector.MCD43C.{}_{}z_{}nm.nc4'.format(outDir,nymd,hour,chd)
    outFilepp   = '{}/aeronet-g5nr.pp.vlidort.vector.MCD43C.{}_{}z_{}nm.nc4'.format(outDir,nymd,hour,chd)
    outFileadd   = '{}/aeronet-g5nr.add.vlidort.vector.MCD43C.{}_{}z.nc4'.format(outDir,nymd,hour)
    
    rcFile   = 'Aod_EOS.rc'
    verbose  = True
    extOnly  = False

    # Initialize VLIDORT class getting aerosol optical properties
    # -----------------------------------------------------------
    vlidort = AERONET_VLIDORT(inFile,rcFile,channel,
                            invFile=invFile,
                            outFileal=outFileal,
                            outFilepp=outFilepp,
                            albedoType=albedoType,
                            brdfFile=brdfFile,
                            ndviFile=ndviFile,
                            lcFile=lcFile,
                            verbose=verbose,
                            extOnly=extOnly,
                            outFileadd=None)

   
    # # Run ext_sampler
    # vlidort.runExt()

    # Run VLIDORT
    # if any(vlidort.nobs):
    #     vlidort.writeNCadd()
        #vlidort.runALMUCANTAR()
        #vlidort.runPRINCIPLE()