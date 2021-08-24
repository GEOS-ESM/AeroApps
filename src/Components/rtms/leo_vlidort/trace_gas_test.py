#!/usr/bin/env python

"""
    Does some simple calculations look at trace gas impacts on RT calculations

    Adapted from benchmark4Amir.py
    Patricia Castellanos, March, 2020

"""

import os
from   netCDF4 import Dataset
import numpy   as np
from MAPL.constants import *
import VLIDORT_POLAR_ 
from polar_vlidort import POLAR_VLIDORT
from benchmark_polar_vlidort import BENCHMARK, MISSING, WrapperFuncs
from scipy.interpolate import interp1d

format   = 'NETCDF4_CLASSIC'
plane_parallel = True


def read_rsr_o3_xsect(inFile):
    f = open(inFile,'r')
    #header
    f.readline()
    f.readline()
    f.readline()
    f.readline()

    wav  = []
    xsec = []
    for l in f:
        w, x = l.split()
        wav.append(w)
        xsec.append(x)

    f.close()

    wav = np.array(wav)
    xsec = np.array(xsec)

    return wav,xsec


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

def get_geom():
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

    return sza,vza,raa,saa,vaa

def get_PTgas_profile(inFile):
    # Read in standard profile
    f = open(inFile,'r')
    f.readline() # header
    f.readline() # header
    pe, te, ze = [],[],[]
    air, o3, o2, h2o, co2, no2 = [], [], [], [], [], []
    for l in f:
        z, p,t,a,O3,O2,H2O,CO2,NO2 = l.split()
        pe.append(float(p))
        te.append(float(t))
        ze.append(float(z))
        air.append(float(a))
        o3.append(float(O3))
        o2.append(float(o2))
        h2o.append(float(H2O))
        co2.append(float(CO2))
        no2.append(float(NO2))

    f.close()
    pe = np.array(pe)
    te = np.array(te)
    ze = np.array(ze)
    air = np.array(air)
    o3  = np.array(o3)
    o2  = np.array(o2)
    h2o = np.array(h2o)
    co2 = np.array(co2)
    no2 = np.array(no2)

    dz = ze[0:-1] - ze[1:]
    zm = zm = ze[1:]+0.5*dz

    # interpoalte to middle
    f = interp1d(ze,air)
    airm = f(zm)
    f = interp1d(ze,o3)
    o3m = f(zm)
    f = interp1d(ze,o2)
    o2m = f(zm)
    f = interp1d(ze,h2o)
    h2om = f(zm)
    f = interp1d(ze,co2)
    co2m = f(zm)
    f = interp1d(ze,no2)
    no2m = f(zm)


    km = pe.shape[0] - 1
    pe.shape = (km+1,1)
    te.shape = (km+1,1)
    ze.shape = (km+1,1)  

    return km,pe,te,ze,airm,o3m,o2m,h2om,co2m,no2m  


def level_1(outFile,Iscan,ixtrack,SZA,VZA,RAA,channels,km,pe,te,ze,U10m,V10m,air,o3,no2):
    #read OCI table
    wav,rod,depol = read_ROD_table()
    frod = interp1d(wav,rod)
    fdepol = interp1d(wav,depol)

    # read xsec
    o3w, o3xsec = read_rsr_o3_xsect('amir/xsect/o3_xsec.txt')

    #albedoType = 'GissCX'
    #albedoType = 'CX'
    albedoType = 'OCIGissCX'
    nstreams = 12

    # AOP vectors
    tau = np.zeros([km,1,1])
    ssa = np.zeros([km,1,1])
    pmom = np.zeros([km,1,1,30,6])

    nscan,nxtrack = SZA.shape
    nscan = 1
    nxtrack = 1
    nch = len(channels)

    # MOP vectors
    gtau = np.zeros([km,1,1])



    # Output vectors to save
    Igranule = np.zeros([nscan,nxtrack,nch])
    Rgranule = np.zeros([nscan,nxtrack,nch])
    BR       = np.zeros([nscan,nxtrack,nch])
    ROD      = np.zeros(nch)

    for ich,channel in enumerate(channels):
        channel = float(channel)
        rcFile   = '../rc/Aod_EOS_{}.rc'.format(int(channel))
        print 'rcFile',rcFile
        verbose  = True

        # Initialize VLIDORT class
        # -----------------------------------------------------------
        vlidort = BENCHMARK(albedoType,channel,outFile,rootDir,rcFile,
                        verbose=verbose,aerosol=False)

        # calculate ROT
        args = [vlidort.channel, pe, ze, te, MISSING, vlidort.verbose]
        vlidortWrapper = WrapperFuncs['ROT_CALC']
        ROT, depol_ratio, rc = vlidortWrapper(*args)
        ROD[ich] = np.sum(ROT)

        # scale ROD to OCI table value
        iwav = channel == wav
        if iwav.sum() == 0:
            #extrapolate for missing wavelengths
            nrod =  frod(channel)
            ndepol = fdepol(channel)
        else:
            nrod = rod[iwav]
            ndepol = depol[iwav]
        ROT = ROT*nrod/ROT.sum()
        depol_ratio = ndepol
        ROD[ich] = np.sum(ROT)

        # get trace gas optical thickness

        # wrapper function based on albedo
        vlidortWrapper = WrapperFuncs[albedoType]

        # run vlidort
        print 'xtrack',ixtrack
        sza = np.array(SZA[Iscan,ixtrack])
        vza = np.array(VZA[Iscan,ixtrack])
        raa = np.array(RAA[Iscan,ixtrack])

        mr   = np.array([1.334])
        args = [channel, nstreams, plane_parallel, ROT, depol_ratio, tau, ssa, pmom,
                pe, ze, te,
                U10m, V10m, mr,
                sza, raa, vza,
                MISSING,
                verbose]

        I, reflectance, surf_reflectance, Q, U, BR_Q, BR_U, rc = vlidortWrapper(*args)

        # save
        Igranule[0,0,ich] = np.squeeze(I)
        Rgranule[0,0,ich] = np.squeeze(reflectance)
        BR[0,0,ich] = np.squeeze(surf_reflectance)
        print 'BR',np.squeeze(surf_reflectance)

    return Igranule,Rgranule,BR,ROD
def writenc(outFile,channels,SZA,SAA,VZA,VAA,I,BR,ROD,Iscan,reflectance=None):
    nscan,nxtrack,nch = I.shape

    # Write data to netcdf file
    nc = Dataset(outFile,'w',format='NETCDF4_CLASSIC')
    nc.title = 'Benchmark case of a Rayleigh plane parallel standard atmosphere'

    dc = nc.createDimension('channels',nch)
    ns = nc.createDimension('number_of_scans',nscan)
    nx = nc.createDimension('ccd_pixels',nxtrack)

    ch = nc.createVariable('channels','f4',('channels',))
    ch.long_name = "Wavelength in nm"
    ch[:] = channels

    rad = nc.createVariable('I','f4',('number_of_scans','ccd_pixels','channels'))
    rad.long_name = "sun normalized TOA radiance"
    rad[:] = I

    if reflectance is not None:
        ref = nc.createVariable('R','f4',('number_of_scans','ccd_pixels','channels'))
        ref.long_name = "TOA reflectance"
        ref[:] = reflectance


    sr = nc.createVariable('surface_reflectance','f4',('number_of_scans','ccd_pixels','channels'))
    sr.long_name = "surface bidirectional reflectance"
    sr[:] = BR  

    a = nc.createVariable('sza','f4',('number_of_scans','ccd_pixels'))
    a.long_name = "solar zenith angle"
    a[:] = SZA[Iscan,:]    

    a = nc.createVariable('saa','f4',('number_of_scans','ccd_pixels'))
    a.long_name = "solar azimiuth angle"
    a[:] = SAA[Iscan,:]        

    a = nc.createVariable('vza','f4',('number_of_scans','ccd_pixels'))
    a.long_name = "sensor zenith angle"
    a[:] = VZA[Iscan,:]    

    a = nc.createVariable('vaa','f4',('number_of_scans','ccd_pixels'))
    a.long_name = "sensor azimiuth angle"
    a[:] = VAA[Iscan,:]        

    rod = nc.createVariable('ROD','f4',('channels',))
    rod.long_name = "Rayleigh Optical Depth"
    rod[:] = ROD

    nc.close()    

#------------------------------------ M A I N ------------------------------------

if __name__ == "__main__":

    outRoot = 'amir/'
    inFile    = '{}/afglms.dat'.format(outRoot)

    # test channels
    #channels = np.array([412.0,443.0,478.0,555.0])
    #nch      = 4
    blue,red,swir = get_channels()
    chdir = {'red': red, 'blue': blue, 'swir': swir}

    rootDir  = '/discover/nobackup/projects/gmao/osse2/pub/c1440_NR/OBS/POLAR_LIDAR/CALIPSO/'
    if not os.path.exists(rootDir):
        rootDir = '/nobackup/3/pcastell/POLAR_LIDAR/CALIPSO/'

    # Read in standard atmosphere
    km, pe, te, ze, air, o3, o2, h2o, co2, no2 =  get_PTgas_profile(inFile)
 
    # read in granule geometry
    SZA,VZA,RAA,SAA,VAA = get_geom()
    Iscan = 600
    Ixtrack = 50

    for band in chdir:
    #for band in ['blue']:
        channels = chdir[band]

        outFile   = '{}/trace_gas_test_{}.nc4'.format(outRoot,band) 
        U10m = np.array([3.0])
        V10m = np.array([4.0])    
        I,reflectance,BR,ROD =  level_1(outFile,Iscan,Ixtrack,SZA,VZA,RAA,channels,km,pe,te,ze,U10m,V10m,air,o3,no2)
        # write to outFile
        writenc(outFile,channels,SZA,SAA,VZA,VAA,I,BR,ROD,Iscan,reflectance=reflectance)    
