#!/usr/bin/env python3

"""
    Does some simple calculations to test trace gas  PACE calculations

    Adapted from benchmark4Amir.py
    Patricia Castellanos, April 2020

"""

import os
import sys
from   netCDF4 import Dataset
from   netCDF4 import Dataset as ncread
import numpy   as np
from MAPL.constants import *
from py_leo_vlidort import VLIDORT_POLAR_
from scipy.interpolate import interp1d
import scipy.integrate as integrate
from pyhdf.SD import SD, SDC
from multiprocessing import Pool

from hyperTest import get_channels,get_geom,get_ROT,get_TOA_unpack,get_TOA
from hyperTest import get_PTWV_profile, writenc

format   = 'NETCDF4_CLASSIC'
MISSING = -1.e+20


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

def get_alpha(A,VMR,rhoe,ze):
    """
    Calculate Absorption optical depth profile
    A - absorption coefficient [m2/molecule]
    VMR - trace gas mixing ratio [vol/vol, dimensionless]
    rhoe - air number density [molecules/m3]
    ze  - profile altitude [m]
    """

    # convert vmr to molecules/m3
    nxe = VMR*rhoe

    # integrate to get the optical depth subcolumns
    km, nbin, nch = A.shape
    alpha = np.zeros([km,nbin,nch])
    for i in range(km):
        for b in range(nbin):
            c1 = A[i,b,:]*nxe[i]
            c2 = A[i,b,:]*nxe[i+1]
            c1.shape = (1, nch)
            c2.shape = (1, nch)
            c  = np.append(c1,c2,axis=0)
            alpha[i,b,:] = np.trapz(c,ze[i:i+2],axis=0)

    #alpha = np.trapz(A*nxe,ze)

    return alpha
#------------------------------------ M A I N ------------------------------------

if __name__ == "__main__":

    outRoot = 'hyperTest/'

    # Pressure [Pa], temperature [K], height [m], water vapor [ppm]  profile - standard atmosphere
    # used to make OCI look up tables
    inFile    = '{}/atrem_tpvmr.nc'.format(outRoot)
    km, pe, te, ze, rhoe, h2oe, DP =  get_PTWV_profile(inFile)

    # Read in Amir's c-k distribution coefficients wavelengths
    # kg has dimensions km,ngbins,nwav.  ngbins = 10
    inFile = 'correlated_k/kg_gas_refined_v2_with_RSR_v0.nc'
    kg_o2, kg_h2o, kg_co, kg_co2, kg_ch4, kg_n2o, g_bins, wav_abs = get_kg(inFile)
    ngbin = len(g_bins)

    # convert kg from pressure space to z space
    # ---------------------------------------------------
    Q = 2.15199993E+25
    C = (Q*28.966) / 6.0225e23 / 1e-6

    # integrate air density in each layer
    rhoint = np.zeros(km)
    for i in range(km):
        rhoint[i] = np.trapz(rhoe[i:i+2],ze[i:i+2])

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
    co_vmr = 0.1*1e-6
    alpha_co = get_alpha(kg_co_z,co_vmr,rhoe,ze)
    
    o2_vmr = 0.21
    alpha_o2 = get_alpha(kg_o2_z,o2_vmr,rhoe,ze)

    co2_vmr = 400.*1.0E-06
    alpha_co2 = get_alpha(kg_co2_z,co2_vmr,rhoe,ze)

    ch4_vmr = 1.8*1.0E-06
    alpha_ch4 = get_alpha(kg_ch4_z,ch4_vmr,rhoe,ze)

    n2o_vmr = 0.3*1.0E-06
    alpha_n2o = get_alpha(kg_n2o_z,n2o_vmr,rhoe,ze)

    sys.exit('stop') 
    # scale water vapor so total precipitable water 
    # is equal to 1 cm
    # ----------------

    # first calculate water vapor column [molecules/m2]
    h2ocol = np.trapz(h2oe*1e-6*rhoe,ze)
    # normalize profile so water vapor column expressed as total precipitable water is equal to 1 cm 
    # use 1 mm of rainfall = 1 kg/m2
    # or 1 cm = 10 kg/m2
    # 10 kg/m2 is equal to 3.34e22 water molecules/cm2
    h2ocolnew = 3.34e22
    # convert to meters
    h2ocolnew = h2ocolnew*1e4
  
    h2oenew = h2oe*(h2ocolnew/h2ocol)
    # get in vmr units
    h2oe_vmr = 1e-6*h2oenew
    alpha_h2o = get_alpha(kg_h2o_z,h2oe_vmr,rhoe,ze)

    # add up all the alphas
    alpha = alpha_h2o + alpha_n2o + alpha_ch4 + alpha_co2 + alpha_o2 + alpha_co

    # append zeros down to wavelength equal to 299.91 nm
    # this is because ozone is not being considered right now
    # eventualy will add O3 cross sections from David Heffner
#    new_wl = np.arange(wav_abs.min()-0.1,299.91,-0.1)
#    all_wl = np.append(wav_abs,new_wl)
    all_wl = wav_abs
#    nnew   = len(new_wl)
#    alpha  = np.append(alpha,np.zeros([km,nnew]),axis=1)

    # flip everything vertically so going from top of atmosphere to surface
    pe = pe[-1::-1]
    te = te[-1::-1]
    ze = ze[-1::-1]
    alpha = alpha[-1::-1,:,:]

    # add dimension to be in km+1,nobs
    pe.shape = (km+1,1)
    te.shape = (km+1,1)
    ze.shape = (km+1,1)


    # read in granule geometry
    Iscan  = 600
    Icross = 1000
    SZA,VZA,RAA,SAA,VAA = get_geom(Iscan,Icross)

    
    # Read in solar irradiance spectrum
    # second dim = wavelength, irradiance
    # units=nm, uW/cm^2/nm
    inFile = '{}/Thuillier_F0.npy'.format(outRoot)
    F0 = np.load(inFile)  
    # interpolate to c-k wavelength bins
    F0_f = interp1d(F0[:,0],F0[:,1],kind='linear',fill_value="extrapolate")
    F0_int = F0_f(all_wl)

    # set up some vlidort stuff
    albedoType = 'OCIGissCX'
    U10m = np.array([3.0])
    V10m = np.array([4.0])
    mr   = 1.334
    nstreams = 12

    # loop through channels
    nproc = 50
    nwl   = len(all_wl)
#    nwl   = 2

    args = []
    ROD  = []
    depol = []
    for ich in np.arange(nwl):
        ch = all_wl[ich]

        # Get Rayleigh
        ROT, depol_ratio = get_ROT(ch,pe,te,ze,verbose=False)           
        ROD.append(np.squeeze(ROT.sum(axis=0)))
        depol.append(depol_ratio)

        for ibin in range(ngbin-1):
            ROT = np.append(ROT,ROT[:,:,0:1],axis=2)
            depol_ratio = np.append(depol_ratio,depol_ratio[0:1])
        
        # trace gas
        alpha_ch = alpha[:,:,ich]
        alpha_ch.shape = (km,1,ngbin)

        # AOP vectors [km,nch,nobs]
        tau = np.zeros([km,ngbin,1])
        ssa = np.zeros([km,ngbin,1])
        pmom = np.zeros([km,ngbin,1,30,6])
 
        # water refractive index
        mr_in = np.ones(ngbin)
        mr_in = mr_in*mr

        # solar irradiance
        F0_in = F0_int[ich]
        F0_in = np.repeat(F0_in,ngbin)
        F0_in.shape = (ngbin,1)

        ch = np.repeat(ch,ngbin)
        args.append([ch,
                     F0_in,
                     ROT,depol_ratio,
                     tau,ssa,pmom,
                     alpha_ch,
                     SZA,VZA,RAA,
                     km,pe,te,ze,
                     nstreams,
                     albedoType,
                     U10m,V10m,mr_in,
                     False])

#       I,reflectance,BR =  get_TOA(ch,
#                                    ROT,depol_ratio,
#                                    tau,ssa,pmom,
#                                    alpha_ch,
#                                    SZA,VZA,RAA,
#                                    km,pe,te,ze,
#                                    nstreams,
#                                    albedoType,
#                                    U10m=U10m,V10m=V10m,mr=np.array([1.334]),
#                                    verbose=False)



    # use multiprocessing
    p = Pool(nproc)
    result = p.map(get_TOA_unpack,args)
    I = []
    reflectance = []
    BR = []
    for r in result:    
        I_r,reflectance_r,BR_r = r
        I.append(np.squeeze(I_r))
        reflectance.append(np.squeeze(reflectance_r))
        BR.append(np.squeeze(BR_r))

    p.close()
    p.join()
    
    # convert to arrays with dimensions nwav,ngbin
    ROD = np.array(ROD)    
    depol_ratio = np.array(depol)

    I  = np.array(I)
    reflectance = np.array(reflectance)
    BR = np.array(BR)
                
    alphaD = alpha.sum(axis=0).T           


    # integrate with g weights to get effective reflectance at each wavelength
    I_eff = np.trapz(I,g_bins,axis=1)
    reflectance_eff = np.trapz(reflectance,g_bins,axis=1)
    BR_eff = np.trapz(BR,g_bins,axis=1)
    alphaD_eff = np.trapz(alphaD,g_bins,axis=1)


    # write to outFile
    outFile   = '{}/outputs/hyperTest_CK_Thuillier_v2.nc4'.format(outRoot)
    writenc(outFile,all_wl,
            F0_int,
            SZA,SAA,VZA,VAA,
            I_eff,BR_eff,reflectance_eff,
            ROD,depol_ratio,
            alphaD_eff,
            pe,te,ze,
            U10m,V10m,mr)    
#

