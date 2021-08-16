#!/usr/bin/env python

"""
    Does some simple calculations to test trace gas  PACE calculations

    Adapted from benchmark4Amir.py
    Patricia Castellanos, April 2020

"""

import os
import sys
from   netCDF4 import Dataset
import numpy   as np
from MAPL.constants import *
from scipy.interpolate import interp1d
import scipy.integrate as integrate
from pyhdf.SD import SD, SDC

format   = 'NETCDF4_CLASSIC'
MISSING = -1.e+20

#---
def get_PT_profile(inFile,model=5):
    """
    Read in height [km],pressure [mb], temperature [K], water vapor vmr profile [ppm]
    for selected model:
    0:Tropical
    1:Mid Latitude Summer
    2:Mid Latitude Winter
    3:Subarctic Summer
    4:Subarctic Winter
    5:US Standard 1962
    6:User Defined Model
    """
    nc = Dataset(inFile)
    # make array because interp doesn't take masked arrays
    pe = np.array(nc.variables['p'][model,:])
    te = np.array(nc.variables['t'][model,:])
    ze = np.array(nc.variables['h'][model,:])
    wv_vmre =np.array( nc.variables['vmr'][model,:])

    km   = len(ze) - 1

    # get DP, from Amir's code
    DP= np.empty([km])
    for i in range(km):
        DP[i] = (pe[i] - pe[i+1])/1013

    # convert mb to Pascal
    pe = pe *100
    DELP = pe[:-1] - pe[1:]
    pm = pe[:-1] - 0.5*DELP

    # get T middles
    f = interp1d(pe,te)
    tm = f(pm)

    # convert km to m
    ze = ze*1000.
    dz = ze[1:] - ze[:-1]
    f = interp1d(pe,ze)
    zm = f(pm)

    # calculate air density 
    g = 9.80616 # gravity m/s2
    R = 8.3145  # gas constant [J mol-1 K-1]
    Na = 6.022e23 # avogadro's number
    Rair = 287.0  # air specific gas constant J kg-1 K-1
    MW_AIR =  28.964*1e-3 # kg/mole 

    AIRDENS = DELP/(dz*g) # kg/m3
    rho = AIRDENS/MW_AIR      # moles/m3
    rho = rho*Na  #[molecules/m3]

    return km, pe, te, ze, rho, wv_vmre, pm, tm, zm, DP, AIRDENS, DELP
#---

def get_RH(pe,te,wv_vmre):
    """
    Convert water vapor vmr [ppm] to relative humidity
    """

    #figure out saturation vapor pressure
    T0 = 273.15  #K
    e0 = 611.    #Pa
    L  = 2.5e6   #J/kg latend heat of vaporization
    Rv = 462.    # J/kg K specific gas constant for water vapor

    tt = (1/T0) - (1/te)
    esat = e0*np.exp(L*tt/Rv)

    # convert VMR to partial pressure
    e = wm_vmre*1e-6*pe

    RH = e/esat

    return RH
#---
def get_kg(inFile):
    """
    Read c-k distribution coefficients from Amir's calculations
    """
    nc = Dataset(inFile)
    # wavelength [nm]
    wav_abs = np.array(nc.variables['wavelengths'][:])
    # c-k  coefficients [not sure about units]
    kg_o2  = nc.variables['kg_o2'][:].T
    kg_h2o = nc.variables['kg_h2o'][:].T
    kg_co  = nc.variables['kg_co'][:].T
    kg_co2 = nc.variables['kg_co2'][:].T
    kg_ch4 = nc.variables['kg_ch4'][:].T
    kg_n2o = nc.variables['kg_n2o'][:].T
    g_bins = nc.variables['g_bins'][:]

    nc.close()

    return kg_o2, kg_h2o, kg_co, kg_co2, kg_ch4, kg_n2o, g_bins, wav_abs

def get_alpha(A,VMR,rho,ze):
    """
    Calculate Absorption optical depth profile
    A - absorption coefficient [m2/molecule]
    VMR - trace gas mixing ratio [vol/vol, dimensionless]
    rho - air number density [molecules/m3]
    ze  - profile altitude [m]
    """

    # convert vmr to molecules/m3
    nx = VMR*rho

    # integrate to get the optical depth subcolumns
    km, nbin, nch = A.shape
    alpha = np.zeros([km,nbin,nch])
    dz = ze[1:]-ze[:-1]
    for i in range(km):
        for b in range(nbin):
            alpha[i,b,:] = A[i,b,:]*nx[i]*dz[i]

    return alpha

def get_abs_o3(inFile):
    """
    get o3 xsection 
    interpolate to ck-bins
    units m2/molecule
    """

    nc = Dataset(inFile)
   
    wav = np.array(wav)
    xsec = np.array(xsec)

    # cross sections are for 1 DU = 2.6967e19 molecules/cm3 @ STP
    # divide by 1DU to get to units cm2/molecule
    xsec = xsec/2.6967e19

    # convert from cm2/molecule to m2/molecule
    xsec = xsec*1e-4

    # append zeros to max ck wavelengths
    wav_new = np.arange(wav.max()+0.1,wav_ck.max()+0.2,0.1)
    wav_o3 = np.append(wav,wav_new)
    nnew = len(wav_new)
    xsec_o3 = np.append(xsec,np.zeros(nnew))

    # interpolate xsec to ck bins wavelengths
    xsec_f = interp1d(wav_o3,xsec_o3,kind='linear')
    xsec_ck = xsec_f(wav_ck)    

    return xsec_ck

#------------------------------------ M A I N ------------------------------------

if __name__ == "__main__":
    atmosFile    = 'hyperTest/atrem_tpvmr.nc'
    ckFile       = 'correlated_k/kg_gas.nc'
    o3File       = 'hyperTest/o3_bremen/xsec_o3_ck_bins_V0_pengwang.nc'
    o3rsrFile    = 'hyperTest/o3_bremen/xsec_o3_rsr_weighted_V0_pengwang.nc'
    solarFile    = 'hyperTest/Thuillier_F0.npy'
    solarFileRSR = 'hyperTest/F0_rsr_weighted_V0.nc'

    co_vmr = 0.1*1e-6
    o2_vmr = 0.21
    co2_vmr = 400.*1.0E-06
    ch4_vmr = 1.8*1.0E-06
    n2o_vmr = 0.3*1.0E-06

    outRoot = 'alphaTable_v0/'
    outFile   = '{}/alpha_CK_Thuillier_o3.nc4'.format(outRoot)

    # make outpath
    if not os.path.exists(outRoot):
        os.makedirs(outRoot)

    # Pressure [Pa], temperature [K], height [m], air density [molecules/m3]  profile - standard atmosphere
    # used to make OCI look up tables
    km, pe, te, ze, rho, wv_vmre, pm, tm, zm, DP, AIRDENS, DELP =  get_PT_profile(atmosFile)

    # Convert H2O VMR profile to RH profile
    #RHe = get_RH(pe,te,wv_vmre)

    # Read in Amir's c-k distribution coefficients wavelengths
    # kg has dimensions km,ngbins,nwav.  ngbins = 10
    kg_o2, kg_h2o, kg_co, kg_co2, kg_ch4, kg_n2o, g_bins, wav_abs = get_kg(ckFile)
    ngbin = len(g_bins)

    # convert kg from pressure space to z space
    # kg - effective absorption coefficient [m2/molecule]
    # ---------------------------------------------------
    Q = 2.15199993E+25
    C = (Q*28.966) / 6.0225e23 / 1e-6

    # integrate air density in each layer
    dz = ze[1:]-ze[:-1]
    rhoint = rho*dz

    DP.shape = (km,1,1)
    rhoint.shape = (km,1,1) 

    for ibin in range(ngbin-1):
        DP = np.append(DP,DP[:,0:1,:],axis=1)
        rhoint = np.append(rhoint,rhoint[:,0:1,:],axis=1)
    kg_o2_z  = kg_o2*C*DP/rhoint
    kg_co_z  = kg_co*C*DP/rhoint
    kg_co2_z = kg_co2*C*DP/rhoint
    kg_ch4_z = kg_ch4*C*DP/rhoint
    kg_n2o_z = kg_n2o*C*DP/rhoint
    kg_h2o_z = kg_h2o*C*DP/rhoint

    # get absorption optical depth with new aborption coefficient
    # dimensions [km,nbin,nch]
    alpha_co = get_alpha(kg_co_z,co_vmr,rho,ze)
    alpha_o2 = get_alpha(kg_o2_z,o2_vmr,rho,ze)
    alpha_co2 = get_alpha(kg_co2_z,co2_vmr,rho,ze)
    alpha_ch4 = get_alpha(kg_ch4_z,ch4_vmr,rho,ze)
    alpha_n2o = get_alpha(kg_n2o_z,n2o_vmr,rho,ze)

    
    # ----
    # OZONE Stuff
    # ---

    # Read in ozone cross sections in m2/molecule
    # interpolated to CK bins
    nc = Dataset(o3File)
    #dimensions [nch,km]
    xsec_o3 = nc.variables['xsec_o3_ck'][:]
#    nc.close()

    # change dimensions to [km,nch]
    xsec_o3 = xsec_o3.T

    # convert from cm2/molecule to m2/molecule
    xsec_o3 = xsec_o3*1e-4  
    nc.close()

    # get rsr weighed for the UV
    nc = Dataset(o3rsrFile)
    wav_rsr = np.array(nc.variables['wav_oci'][:])
    xsec_o3_rsr = nc.variables['xsec_o3_rsr'][:]

    # chanege dimesions to [km,ch]
    xsec_o3_rsr = xsec_o3_rsr.T

    # convert from cm2/molecule to m2/molecule
    xsec_o3_rsr = xsec_o3_rsr*1e-4

    # Read in solar irradiance spectrum
    # dims = [wavelength,2]
    # second dim = wavelength(0), irradiance(1)
    # units=nm, uW/cm^2/nm
    F0 = np.load(solarFile)
    # interpolate to c-k wavelength bins
    F0_f = interp1d(F0[:,0],F0[:,1],kind='linear',fill_value="extrapolate")
    F0_int = F0_f(wav_abs)

    # Convert to W/m2/nm
    F0_int = F0_int*1e-6*1e4

    # append together UV-rsr weight with C-K interpolated
    # only OCI channels not covered by CK-bins
    # for v0 RSR, ck min = 575, i=53
    # for 2.5 nm RSR, ck min = 570, i=103
    # wav_abs is not strictly monotinic anymore
    iend = 53
    wav_abs = np.append(wav_rsr[:iend],wav_abs)

    xsec_o3 = np.append(xsec_o3_rsr[:,:iend],xsec_o3,axis=1)

    # append RSR weighted solar irradiance
    nc = Dataset(solarFileRSR)
    F0_rsr = np.array(nc.variables['F0_rsr'][:])
    # Convert to W/m2/nm
    F0_rsr = F0_rsr*1e-6*1e4

    F0_int = np.append(F0_rsr[:iend],F0_int)

    # set up an nbins variable.  
    # C-K wavelengths have 10
    # UV have 1
    nbins = np.ones(len(wav_abs))
    nbins[53:] = ngbin

    # add zeros to beginning of alpha
    aa = np.zeros([km,ngbin,iend])
    alpha_n2o = np.append(aa,alpha_n2o,axis=2)
    alpha_ch4 = np.append(aa,alpha_ch4,axis=2)
    alpha_co2 = np.append(aa,alpha_co2,axis=2)
    alpha_o2  = np.append(aa,alpha_o2,axis=2)
    alpha_co  = np.append(aa,alpha_co,axis=2)

    # append zeros to begining of h2o cross section
    kg_h2o_z = np.append(aa,kg_h2o_z,axis=2)
    kg_co_z = np.append(aa,kg_co_z,axis=2)
    kg_co2_z = np.append(aa,kg_co2_z,axis=2)

    # add up all the alphas
    alpha = alpha_n2o + alpha_ch4 + alpha_co2 + alpha_o2 + alpha_co
    
    # flip everything vertically so going from top of atmosphere to surface
    pe    = pe[-1::-1]
    te    = te[-1::-1]
    ze    = ze[-1::-1]
    pm    = pm[-1::-1]
    tm    = tm[-1::-1]
    zm    = zm[-1::-1]
    AIRDENS = AIRDENS[-1::-1]
    DELP    = DELP[-1::-1]
    rho  = rho[-1::-1]
    alpha = alpha[-1::-1,:,:]
    alpha_n2o = alpha_n2o[-1::-1,:,:]
    alpha_ch4 = alpha_ch4[-1::-1,:,:]
    alpha_co2 = alpha_co2[-1::-1,:,:]
    alpha_o2 = alpha_o2[-1::-1,:,:]
    alpha_co = alpha_co[-1::-1,:,:]
    xsec_o3 = xsec_o3[-1::-1,:]
    kg_h2o_z = kg_h2o_z[-1::-1,:,:]

    # write to outFile
    nch = len(wav_abs)

    nc = Dataset(outFile,'w',format='NETCDF4_CLASSIC')
    nc.title = 'absorption optical depths and cross sections for C-K wavelengths'

    dc = nc.createDimension('channels',nch)
    dck = nc.createDimension('ckbin',ngbin)
    dk = nc.createDimension('lev',km)
    de = nc.createDimension('leve',km+1)

    ch = nc.createVariable('channels','f4',('channels',))
    ch.long_name = "Wavelength in nm"
    ch[:] = wav_abs

    b  = nc.createVariable('nbins','i4',('channels',))
    b.long_name = 'number of sub-bins for each channel'
    b[:] = nbins

    g  = nc.createVariable('g_bins','f4',('ckbin',))
    g.long_name = 'C-K g integration factor'
    g[:] = g_bins

    f  = nc.createVariable('solar_irradiance','f4',('channels',))
    f.long_name = 'Thuillier solar irradiance spectrum'
    f.units = 'uW/cm^2/nm'
    f[:] = F0_int 

    ke = nc.createVariable('PE','f4',('leve',))
    ke.long_name = 'Pressure at layer edge'
    ke.units     = 'Pa'
    ke[:] = pe

    ke = nc.createVariable('TE','f4',('leve',))
    ke.long_name = 'Temperature at layer edge'
    ke.units     = 'K'
    ke[:] = te

    ke = nc.createVariable('ZE','f4',('leve',))
    ke.long_name = 'height above surface'
    ke.units = 'm'
    ke[:] = ze

    ke = nc.createVariable('P','f4',('lev',))
    ke.long_name = 'Pressure at layer middle'
    ke.units     = 'Pa'
    ke[:] = pm

    ke = nc.createVariable('T','f4',('lev',))
    ke.long_name = 'Temperature at layer middle'
    ke.units     = 'K'
    ke[:] = tm

    ke = nc.createVariable('Z','f4',('lev',))
    ke.long_name = 'height above surface at layer middle'
    ke.units = 'm'
    ke[:] = zm

    ke = nc.createVariable('AIRDENS','f4',('lev',))
    ke.long_name = 'air density'
    ke.units = 'kg/m3'
    ke[:] = AIRDENS

    ke = nc.createVariable('DELP','f4',('lev',))
    ke.long_name = 'layer pressure thickness'
    ke.units = 'Pa'
    ke[:] = DELP

    ke = nc.createVariable('xsec_o3','f4',('lev','channels',))
    ke.long_name = 'ozone cross section profile'
    ke.units = 'm2/molecule'
    ke[:] = xsec_o3

    ke = nc.createVariable('xsec_h2o','f4',('lev','ckbin','channels',))
    ke.long_name = 'water vapor cross section profile'
    ke.units = 'm2/molecule'
    ke[:] = kg_h2o_z

    ke = nc.createVariable('xsec_co','f4',('lev','ckbin','channels',))
    ke.long_name = 'CO cross section profile'
    ke.units = 'm2/molecule'
    ke[:] = kg_co_z

    ke = nc.createVariable('xsec_co2','f4',('lev','ckbin','channels',))
    ke.long_name = 'CO2 cross section profile'
    ke.units = 'm2/molecule'
    ke[:] = kg_co2_z


    ke = nc.createVariable('alpha_total','f4',('lev','ckbin','channels',))
    ke.long_name = 'total trace gas absorption optical depth (excludes water vapor & ozone)'
    ke.units = 'None'
    ke[:] = alpha

    ke = nc.createVariable('alpha_co','f4',('lev','ckbin','channels',))
    ke.long_name = 'CO absorption optical depth'
    ke.units = 'None'
    ke[:] = alpha_co

    ke = nc.createVariable('alpha_o2','f4',('lev','ckbin','channels',))
    ke.long_name = 'O2 absorption optical depth'
    ke.units = 'None'
    ke[:] = alpha_o2

    ke = nc.createVariable('alpha_co2','f4',('lev','ckbin','channels',))
    ke.long_name = 'CO2 absorption optical depth'
    ke.units = 'None'
    ke[:] = alpha_co2

    ke = nc.createVariable('alpha_ch4','f4',('lev','ckbin','channels',))
    ke.long_name = 'CH4 absorption optical depth'
    ke.units = 'None'
    ke[:] = alpha_ch4

    ke = nc.createVariable('alpha_n2o','f4',('lev','ckbin','channels',))
    ke.long_name = 'N2O absorption optical depth'
    ke.units = 'None'
    ke[:] = alpha_n2o

    nc.close()
