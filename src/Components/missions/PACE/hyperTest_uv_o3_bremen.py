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

from hyperTest_o3 import get_channels,get_geom,get_ROT,get_TOA_unpack,get_TOA
from hyperTest_o3 import get_PTWV_profile, get_o3_profile, writenc

format   = 'NETCDF4_CLASSIC'
MISSING = -1.e+20


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

def get_abs_o3_rsr(inFile,te):
    """
    get o3 xsection
    use values that are RSR weighted
    interpolate to temp profile te
    units m2/molecule
    """
    # rsr weighted xsec
    nc = Dataset(inFile)
    wav_rsr = np.array(nc.variables['wav_oci'][:])
    xsec_rsr = nc.variables['xsec_o3_rsr'][:]
    temp    = nc.variables['temp'][:]
    nc.close()

    # get xsecstion for te temperature profile
    xsec_f = interp1d(temp,xsec_rsr,kind='linear',fill_value='extrapolate')

    # dimensions ke,nwav
    xsec_int = xsec_f(te).T

    # convert from cm2/molecule to m2/molecule
    xsec_int = xsec_int*1e-4

    return wav_rsr,xsec_int

def get_rod_rsr(inFile,ch):
    nc = Dataset(inFile)
    rod = nc.variables['rod_rsr'][:]
    depol = nc.variables['depol_rsr'][:]
    wav = nc.variables['wav_oci'][:]
    nc.close()

    i = wav == ch

    return rod[i],depol[i]
#------------------------------------ M A I N ------------------------------------

if __name__ == "__main__":

    outRoot = 'hyperTest/'
    outFile   = '{}/outputs/hyperTest_UV_Thuillier_o3.nc4'.format(outRoot)

    # Pressure [Pa], temperature [K], height [m], water vapor [ppm]  profile - standard atmosphere
    # used to make OCI look up tables
    inFile    = '{}/atrem_tpvmr.nc'.format(outRoot)
    km, pe, te, ze, rhoe, h2oe, DP =  get_PTWV_profile(inFile)

    # ----
    # OZONE Stuff
    # ---
    # Read in mcclams ozone mixing ratio profile [unitless]
    inFile = '{}/mcclams.dat'.format(outRoot)
    z_o3, o3_vmr = get_o3_profile(inFile)

    # interpolate o3_vmr to PT profile
    o3_f = interp1d(z_o3,o3_vmr,kind='linear',fill_value="extrapolate")
    o3_vmre = o3_f(ze)

    # Read in ozone cross sections
    # weighted by rsr
    # also interpolated to te profile
    inFile = '{}/o3_bremen/xsec_o3_rsr_weighted_V0.nc'.format(outRoot)
    wav_abs, abs_o3_z = get_abs_o3_rsr(inFile,te)
    
    # integrate conc_o3*xsec_o3 to get o3  alpha
    nwav = len(wav_abs)
    alpha_o3 = np.zeros([km,nwav])
    for i in range(km):
        o3_conc = o3_vmre[i:i+2]*rhoe[i:i+2]
        o3_conc.shape = (2,1)
        alpha_o3[i,:] = np.trapz(o3_conc*abs_o3_z[i:i+2,:],ze[i:i+2],axis=0)
    
    alpha = alpha_o3
    all_wl = wav_abs

    # flip everything vertically so going from top of atmosphere to surface
    pe = pe[-1::-1]
    te = te[-1::-1]
    ze = ze[-1::-1]
    alpha = alpha[-1::-1,:]

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
    nproc = 20
    nwl   = len(all_wl)
#    nwl   = 2

    args = []
    ROD  = []
    depol = []
    nch = 1
    for ich in np.arange(nwl):
        ch = all_wl[ich]

        # Get Rayleigh
        ROT, depol_ratio = get_ROT(ch,pe,te,ze,verbose=False)

        # Read in RSR weighted ROD and depol
        inFile = '{}/rod_rsr_weighted_V0.nc'.format(outRoot)
        rod,depol_ratio = get_rod_rsr(inFile,ch)

        # Scale ROT to ROD
        ROT = ROT*rod/ROT.sum()

        ROD.append(rod)
        depol.append(depol_ratio)

        # trace gas
        alpha_ch = alpha[:,ich]
        alpha_ch.shape = (km,1,nch)

        # AOP vectors [km,nch,nobs]
        tau = np.zeros([km,nch,1])
        ssa = np.zeros([km,nch,1])
        pmom = np.zeros([km,nch,1,30,6])
 
        # water refractive index
        mr_in = np.ones(nch)
        mr_in = mr_in*mr

        # solar irradiance
        F0_in = F0_int[ich:ich+1]
        F0_in.shape = (nch,1)

        
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
    
    # concatenate arrays
    ROD = np.concatenate(ROD)
    depol_ratio = np.concatenate(depol)

    I  = np.array(I)
    reflectance = np.array(reflectance)
    BR = np.array(BR)

    alphaD = alpha.sum(axis=0)
    # write to outFile
    writenc(outFile,all_wl,
            F0_int,
            SZA,SAA,VZA,VAA,
            I,BR,reflectance,
            ROD,depol_ratio,
            alphaD,
            pe,te,ze,
            U10m,V10m,mr)    
#

