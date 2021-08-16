#!/usr/bin/env python
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

#------------------------------------ M A I N ------------------------------------

if __name__ == "__main__":
    # Read in solar irradiance spectrum
    # second dim = wavelength, irradiance
    # units=nm, uW/cm^2/nm
    inFile = 'Thuillier_F0.npy'
    F0 = np.load(inFile)

    # read in RSR
    #rsrFile = 'OCI_RSR_v0.hdf'
    rsrFile = 'OCI_RSR_2.5nm_v7_amended.hdf'
    rsr, wav_rsr, wav_oci = get_rsr(rsrFile)
    noci = len(wav_oci)

    # interpolate rsr to F0 wavelengths
    rsr_f = interp1d(wav_rsr,rsr,kind='linear',fill_value="extrapolate")
    rsr_int = rsr_f(F0[:,0])

    # smooth F0
    f0_smooth = np.zeros([noci])
    for ich in range(noci):
        norm = np.trapz(rsr_int[ich,:],F0[:,0])
        if norm != 0:
            f0_smooth[ich] = np.trapz(F0[:,1]*rsr_int[ich,:],F0[:,0])/norm

    # write to a netcdf file
    #outFile = 'F0_rsr_weighted_V0.nc'
    outFile = 'F0_rsr_weighted_V7.nc'
    nc = Dataset(outFile,'w')
    nc.F0File = 'Thuillier_F0.npy'
    nc.rsrFile = rsrFile
    nc.contact = 'Patricia Castellanos <patricia.castellanos@nasa.gov>'
    nc.comments = 'solar irradiance weighted with OCI RSR'

    dlbl = nc.createDimension('wav_f0',len(F0[:,0]))
    doci = nc.createDimension('wav_oci',noci)

    ch = nc.createVariable('wav_f0','f4',('wav_f0',))
    ch.long_name = 'solar irradiance wavelength'
    ch.units = 'nm'
    ch[:] = F0[:,0]

    ch = nc.createVariable('wav_oci','f4',('wav_oci',))
    ch.long_name = 'oci wavelengths'
    ch.units = 'nm'
    ch[:] = wav_oci

    x = nc.createVariable('F0','f4',('wav_f0',))
    x.long_name = 'solar irradiance'
    x.units = 'uW/cm^2/nm'
    x[:] = F0[:,1]

    x = nc.createVariable('F0_rsr','f4',('wav_oci',))
    x.long_name = 'solar irradiance weighted by OCI RSR'
    x.units = 'uW/cm^2/nm'
    x[:] = f0_smooth

    nc.close()
