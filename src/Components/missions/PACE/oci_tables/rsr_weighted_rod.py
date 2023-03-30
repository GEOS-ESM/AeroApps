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

# --
def read_ROD_table(inFile):
    f = open(inFile)

    for i in range(16):
        f.readline() #header

    wav = []
    rod = []
    depol = []
    for l in f:
        w, r, d = l.split()
        wav.append(w)
        rod.append(r)
        depol.append(d)

    f.close()

    wav = np.array(wav).astype('float')
    rod = np.array(rod).astype('float')
    depol = np.array(depol).astype('float')

    return wav, rod, depol

#------------------------------------ M A I N ------------------------------------

if __name__ == "__main__":
    # read xsec
    inFile = 'rayleigh_bodhaine.txt'
    wav_lbl,rod,depol = read_ROD_table(inFile)

    # read in RSR
    #rsrFile = '../hyperTest/OCI_RSR_v0.hdf'
    rsrFile = '../hyperTest/OCI_RSR_2.5nm_v7_amended.hdf'
    rsr, wav_rsr, wav_oci = get_rsr(rsrFile)
    noci = len(wav_oci)

    # interpolate rsr to LBL wavelengths
    rsr_f = interp1d(wav_rsr,rsr,kind='linear',fill_value="extrapolate")
    rsr_int = rsr_f(wav_lbl)

    # smooth ROD
    rod_smooth = np.zeros([noci])
    depol_smooth = np.zeros([noci])
    for ich in range(noci):
        norm = np.trapz(rsr_int[ich,:],wav_lbl)
        if norm != 0:
            rod_smooth[ich] = np.trapz(rod*rsr_int[ich,:],wav_lbl)/norm
            depol_smooth[ich] = np.trapz(depol*rsr_int[ich,:],wav_lbl)/norm

    # write to a netcdf file
    #outFile = 'rod_rsr_weighted_V0.nc'
    outFile = 'rod_rsr_weighted_V7.nc'
    nc = Dataset(outFile,'w')
    nc.rodFile = 'rayleigh_bodhaine.txt'
    nc.rsrFile = rsrFile
    nc.contact = 'Patricia Castellanos <patricia.castellanos@nasa.gov>'
    nc.comments = 'rayleigh optical depth and depol ration weighted with OCI RSR'

    dlbl = nc.createDimension('wav_lbl',len(wav_lbl))
    doci = nc.createDimension('wav_oci',noci)

    ch = nc.createVariable('wav_lbl','f4',('wav_lbl',))
    ch.long_name = 'lbl wavelength'
    ch.units = 'nm'
    ch[:] = wav_lbl

    ch = nc.createVariable('wav_oci','f4',('wav_oci',))
    ch.long_name = 'oci wavelengths'
    ch.units = 'nm'
    ch[:] = wav_oci

    x = nc.createVariable('rod','f4',('wav_lbl',))
    x.long_name = 'lbl ROD'
    x.units = 'unitless'
    x[:] = rod

    x = nc.createVariable('depol','f4',('wav_lbl',))
    x.long_name = 'lbl depolarizaton ratio'
    x.units = 'unitless'
    x[:] = depol

    x = nc.createVariable('rod_rsr','f4',('wav_oci',))
    x.long_name = 'ROD weighted by OCI RSR'
    x.units = 'unitless'
    x[:] = rod_smooth

    x = nc.createVariable('depol_rsr','f4',('wav_oci',))
    x.long_name = 'depolarization ratio weighted by OCI RSR'
    x.units = 'unitless'
    x[:] = depol_smooth

    nc.close()
