#!/usr/bin/env python

"""
    Does some simple calculations to benchmark PACE calculations

    Adapted from polar_vlidort
    Patricia Castellanos, May, 2017

"""

import os
from   netCDF4 import Dataset
import numpy   as np
from MAPL.constants import *
import VLIDORT_POLAR_ 
from polar_vlidort import POLAR_VLIDORT
from benchmark_polar_vlidort import BENCHMARK, MISSING, WrapperFuncs
from scipy.interpolate import interp1d
from   mieobs  import getAOPvector, getEdgeVars, getAOPint, getAOPext

format   = 'NETCDF4_CLASSIC'
plane_parallel = True


# Generic Lists of Varnames and Units
VNAMES_DU = ['DU001','DU002','DU003','DU004','DU005']
VNAMES_SS = ['SS001','SS002','SS003','SS004','SS005']
VNAMES_BC = ['BCPHOBIC','BCPHILIC']
VNAMES_OC = ['OCPHOBIC','OCPHILIC']
VNAMES_SU = ['SO4']
MieVarsNames = ['ext','scatext','backscat','aback_sfc','aback_toa','depol','ext2back','tau','ssa','g']
MieVarsUnits = ['km-1','km-1','km-1 sr-1','sr-1','sr-1','unitless','sr','unitless','unitless','unitless']

META    = ['DELP','PS','RH','AIRDENS']
AERNAMES = VNAMES_SU + VNAMES_SS + VNAMES_OC + VNAMES_BC + VNAMES_DU
SDS_AER = META + AERNAMES
MISSING = -1.e+20

i = 1547
j = 1032

sleave = np.array([0.0078952684998512])
sleave.shape = (1,) + sleave.shape

class GEOS(object):
    def __init__(self,inFile,metFile):
        self.inFile = inFile
        self.metFile = metFile
        self.SDS_AER = SDS_AER

        self.readSampledGEOS()

        pe, ze, te = getEdgeVars(self)

        self.pe = pe # (km,nobs)
        self.ze = ze 
        self.te = te

        km,nobs = pe.shape
        self.km = km



    #---
    def readSampledGEOS(self):
        """
        Read in model sampled track
        """
        nc       = Dataset(self.inFile)
        SDS = self.SDS_AER
        for sds in SDS:
            self.__dict__[sds] = []
            if len(nc.variables[sds][:].shape) == 3:
                var = np.array([nc.variables[sds][0,i,j]])
            else:
                var = nc.variables[sds][0,:,i,j]

            var.shape = (1,) + var.shape
            self.__dict__[sds] = var   

        nc.close()

        nc = Dataset(self.metFile) 
        for sds in ['U10M','V10M']:
            var = np.array([nc.variables[sds][0,i,j]])
            self.__dict__[sds] = var   


def read_ROD_table():
    f = open('amir/OCI_ROD_Table_adjusted.txt','r')
    f.readline() #header
    wav = f.readline().split()
    f.readline() #wav center
    f.readline() #wav width
    f.readline() #F0
    rod = f.readline().split()
    depol = f.readline().split() 

    f.close()

    rod = rod[2:]
    rod = np.array(rod).astype('float')

    depol = depol[2:]
    depol = np.array(depol).astype('float')

    wav = wav[3:]
    wav = np.array(wav).astype('float')

    return wav,rod,depol


def get_channels():
    inFile = '/nobackup/PACE/L1B/Y2020/M03/D24/OCI2020084005000.L1B_PACE.nc'
    nc = Dataset(inFile)
    grp = nc.groups['sensor_band_parameters']
    blue = grp.variables['blue_wavelength'][:]
    red  = grp.variables['red_wavelength'][:]
    swir = grp.variables['SWIR_wavelength'][:]

    return blue,red,swir

def get_geom(i,j):
    inFile = '/nobackup/PACE/L1B/Y2020/M03/D24/OCI2020084005000.L1B_PACE.nc'
    nc = Dataset(inFile)
    grp = nc.groups['geolocation_data']
    sza = grp.variables['solar_zenith'][:]
    saa = grp.variables['solar_azimuth'][:]
    vza = grp.variables['sensor_zenith'][:]
    vaa = grp.variables['sensor_azimuth'][:]

    # make azimuths clockwise from north
    I = saa < 0
    saa[I] = 360. + saa[I]

    I = vaa < 0
    vaa[I] = 360. + vaa[I]

    # define RAA according to photon travel direction
    saa = saa + 180.0
    I = saa >= 360.
    saa[I] = saa[I] - 360.

    raa = vaa - saa   

    I = raa < 0
    raa[I] = raa[I] + 360.0

    saa = grp.variables['solar_azimuth'][:]
    vaa = grp.variables['sensor_azimuth'][:]

    return sza[i,j],vza[i,j],raa[i,j],saa[i,j],vaa[i,j]

def level_1(geos,SZA,VZA,RAA,channel,rcFile):
    #read OCI table
    wav,rod,depol = read_ROD_table()
    frod = interp1d(wav,rod)
    fdepol = interp1d(wav,depol)

    albedoType = 'OCIGissCX_NOBM_CLOUD'
    nstreams = 12

    # AOP vectors
    tau,ssa,g,pmom = getAOPvector(geos,channel,
                         vnames=AERNAMES,
                         Verbose=True,
                         rcfile=rcFile,
                         nMom=300)

    tauI = np.zeros(tau.shape)
    ssaI = np.zeros(ssa.shape)
    pmomI = np.zeros(pmom.shape)

    tauL = tauI
    ssaL = ssaI
    pmomL = pmomI

    # Initialize VLIDORT class 
    # -----------------------------------------------------------
    outFile = 'None'
    rootDir = '/nobackup/3/pcastell/POLAR_LIDAR/CALIPSO/'
    verbose = True
    vlidort = BENCHMARK(albedoType,channel,outFile,rootDir,rcFile,
                    verbose=verbose,aerosol=False)

    args = [channel, geos.pe, geos.ze, geos.te, MISSING, verbose]
    vlidortWrapper = WrapperFuncs['ROT_CALC']
    ROT, depol_ratio, rc = vlidortWrapper(*args) 
    ROD = np.sum(ROT)

    print 'starting ROD',ROD

    # wrapper function based on albedo
    vlidortWrapper = WrapperFuncs[albedoType]

    # scale to OCI table value
    iwav = channel == wav
    if iwav.sum() == 0:
        #extrapolate for missing wavelengths
        nrod =  frod(channel)
        ndepol = fdepol(channel)
    else:
        nrod = rod[iwav]
        ndepol = depol[iwav]
    ROT = ROT*nrod/ROT.sum()
    print 'standard rod',nrod
    # use surface pressure for scaling
    ROT = ROT*geos.pe[-1,0]/101300.0
    depol_ratio = ndepol
    ROD = np.sum(ROT)



    mr   = np.array([1.334])
    args = [channel, nstreams, plane_parallel, ROT, depol_ratio, tau, ssa, pmom,
            tauI, ssaI, pmomI, tauL, ssaL, pmomL,
            geos.pe, geos.ze, geos.te,
            geos.U10M, geos.V10M, mr, sleave,
            SZA, RAA, VZA,
            MISSING,
            verbose]

    I, reflectance, surf_reflectance, Q, U, BR_Q, BR_U, rc = vlidortWrapper(*args)    

    geos.rot = ROT
    geos.ssa = ssa
    geos.tau = tau
    return I,reflectance,surf_reflectance,ROD



#------------------------------------ M A I N ------------------------------------

if __name__ == "__main__":

    aerFile   = '/nobackup/PACE/LevelB/Y2006/M03/D24/pace-g5nr.lb.aer_Nv.20060324_005000.nc4'
    metFile   = '/nobackup/PACE/LevelB/Y2006/M03/D24/pace-g5nr.lb.met_Nv.20060324_005000.nc4'
    channel   = 445
    geos = GEOS(aerFile,metFile)
    # read in granule geometry
    SZA,VZA,RAA,SAA,VAA = get_geom(i,j)


    U10m = np.array([3.0])
    V10m = np.array([4.0])    
    rcFile = 'rc/Aod_EOS.{}.rc'.format(channel)
    I,reflectance,BR,ROD =  level_1(geos,SZA,VZA,RAA,channel,rcFile)
