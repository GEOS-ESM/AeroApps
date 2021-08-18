#!/usr/bin/env python
"""
interpolate o3 spectra to C-K bin centers
"""

import os
import sys
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import numpy as np
from netCDF4 import Dataset


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

    # read in c-k wavelengths [nm]
    inFile = '../../correlated_k/kg_gas.nc'
    nc = Dataset(inFile)
    wav_ck = np.array(nc.variables['wavelengths'][:])
    nck    = len(wav_ck)

    # interpolate ozone xsec to ck wavelengths
    # when you go past the last wavelength, use zeros
    xsec_ck = np.zeros([nck,km])
    for itemp in range(km):
        xsec_f = interp1d(wav_o3,xsec_o3[:,itemp],kind='linear',fill_value=0.0,bounds_error=False)
        xsec_ck[:,itemp] = xsec_f(wav_ck)

    # write to a netcdf file
    outFile = 'xsec_o3_ck_bins_V0_pengwang.nc'
    nc = Dataset(outFile,'w')
    nc.xsecFile = 'Ozone_abs_x_wTemperatureFit.dat'
    nc.tempFile = 'atrem_tpvmr.nc'
    nc.ckFile = 'kg_gas.nc'
    nc.contact = 'Patricia Castellanos <patricia.castellanos@nasa.gov>'
    nc.comments = 'ozone cross section profiles interpolated to correlated-k wavelength bins'

    do3 = nc.createDimension('wav_o3',len(wav_o3))
    doci = nc.createDimension('wav_ck',nck)
    dtemp = nc.createDimension('lev',km)

    ch = nc.createVariable('wav_o3','f4',('wav_o3',))
    ch.long_name = 'o3 xsec wavelength (high resolution)'
    ch.units = 'nm'
    ch[:] = wav_o3

    ch = nc.createVariable('wav_ck','f4',('wav_ck',))
    ch.long_name = 'correlated-k wavelengths'
    ch.units = 'nm'
    ch[:] = wav_ck

    t = nc.createVariable('temp','f4',('lev',))
    t.long_name = 'temperature'
    t.units = 'K'
    t[:] = tm

    x = nc.createVariable('xsec_o3','f4',('wav_o3','lev',))
    x.long_name = 'o3 cross section (original)'
    x.units = 'cm2/molecule'
    x[:] = xsec_o3

    x = nc.createVariable('xsec_o3_ck','f4',('wav_ck','lev'))
    x.long_name = 'o3 cross section interpolated to ck bins'
    x.units = 'cm2/molecule'
    x[:] = xsec_ck

    nc.close()
