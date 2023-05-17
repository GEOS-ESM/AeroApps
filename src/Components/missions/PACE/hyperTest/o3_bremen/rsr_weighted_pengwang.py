#!/usr/bin/env python3
"""
read and plot o3 spectra from text file
"""

import os
import sys
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import numpy as np
from netCDF4 import Dataset
from pyhdf.SD import SD, SDC

def read_o3(inFile,te):
    f = open(inFile)
    nhead = 10
    for i in range(nhead):
        hh = f.readline()

    wav_o3, c0, c1, c2 = [],[],[],[]
    for l in f:
        a = np.array(l.split()).astype(float)
        wav_o3.append(a[0])
        c0.append(a[1])
        c1.append(a[2])
        c2.append(a[3])
    f.close()

    wav_o3 = np.array(wav_o3)
    c0     = np.array(c0)
    c1     = np.array(c1)
    c2     = np.array(c2)
    T0     = 273.15

    # calculate xsec for te
    nwav = len(wav_o3)
    km   = len(te)
    xsec_o3 = np.zeros([nwav,km])

    for i,t in enumerate(te):
        xsec_o3[:,i] = c0 + c1*(t-T0) + c2*(t-T0)**2

    xsec_o3 = xsec_o3*1e-20

    return wav_o3, xsec_o3

#---
def get_rsr(inFile):
    """
    Read in OCI RSR File
    """
    hdf = SD(inFile, SDC.READ)
    rsr = hdf.select('RSR')[:]
    wav_rsr = hdf.select('rsrwave')[:]
    wav_oci = hdf.select('wave')[:]
    hdf.end()

    return rsr.T, wav_rsr, wav_oci

#---
def get_T_profile(inFile,model=5):
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
    te = np.array(nc.variables['t'][model,:])
    ze = np.array(nc.variables['h'][model,:])
    pe = np.array(nc.variables['p'][model,:])

    km = len(te) - 1

    # get T middles
    DELP = pe[:-1] - pe[1:]
    pm = pe[:-1] - 0.5*DELP
    f = interp1d(pe,te)
    tm = f(pm)

    return km,tm


#------------------------------------ M A I N ------------------------------------

if __name__ == "__main__":
    # read T profile
    inFile = '../atrem_tpvmr.nc'
    km,tm = get_T_profile(inFile)
    
    # read xsec
    inFile = 'Ozone_abs_x_wTemperatureFit.dat'
    # wav [nm], C0, C1(T), C2(T^2)
    # xsec is in cm2/molecule
    wav_o3,xsec_o3 = read_o3(inFile,tm)

    # read in RSR
    inFile = '../OCI_RSR_v0.hdf'
    rsr, wav_rsr, wav_oci = get_rsr(inFile)
    noci = len(wav_oci)

    # extent xsec_o3 to max rsr wavelength with zeros
    wav_new = np.arange(wav_o3.max()+0.1,wav_rsr.max()+0.1,0.1)
    nnew = len(wav_new)
    wav_o3 = np.append(wav_o3,wav_new)
    xsec_o3 = np.append(xsec_o3,np.zeros((nnew,km)),axis=0)

    # interpolate rsr to ozone wavelengths
    rsr_f = interp1d(wav_rsr,rsr,kind='linear',fill_value="extrapolate")
    rsr_int = rsr_f(wav_o3)

    # smooth xsection
    xsec_smooth = np.zeros([noci,km])
    for ich in range(noci):
        norm = np.trapz(rsr_int[ich,:],wav_o3)
        if norm != 0:
            for itemp in range(km):
                xsec_smooth[ich,itemp] = np.trapz(xsec_o3[:,itemp]*rsr_int[ich,:],wav_o3)/norm

    # write to a netcdf file
    outFile = 'xsec_o3_rsr_weighted_V0_pengwang.nc'
    nc = Dataset(outFile,'w')
    nc.xsecFile = 'Ozone_abs_x_wTemperatureFit.dat'
    nc.rsrFile = 'OCI_RSR_v0.hdf'
    nc.tempFile = 'atrem_tpvmr.nc'
    nc.contact = 'Patricia Castellanos <patricia.castellanos@nasa.gov>'
    nc.comments = 'ozone cross section profiles weighted with OCI RSR'

    do3 = nc.createDimension('wav_o3',len(wav_o3))
    doci = nc.createDimension('wav_oci',noci)
    dtemp = nc.createDimension('lev',km)

    ch = nc.createVariable('wav_o3','f4',('wav_o3',))
    ch.long_name = 'o3 xsec wavelength (original)'
    ch.units = 'nm'
    ch[:] = wav_o3

    ch = nc.createVariable('wav_oci','f4',('wav_oci',))
    ch.long_name = 'oci wavelengths'
    ch.units = 'nm'
    ch[:] = wav_oci

    t = nc.createVariable('temp','f4',('lev',))
    t.long_name = 'temperature profile'
    t.units = 'K'
    t[:] = tm

    x = nc.createVariable('xsec_o3','f4',('wav_o3','lev',))
    x.long_name = 'o3 cross section (high-res)'
    x.units = 'cm2/molecule'
    x[:] = xsec_o3

    x = nc.createVariable('xsec_o3_rsr','f4',('wav_oci','lev'))
    x.long_name = 'o3 cross section weighted by OCI RSR'
    x.units = 'cm2/molecule'
    x[:] = xsec_smooth

    nc.close()
