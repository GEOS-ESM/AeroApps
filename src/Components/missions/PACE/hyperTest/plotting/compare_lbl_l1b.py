#!/usr/bin/env python

"""
Compare LBL and C-K gas transmittance
"""

import os
import sys
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from pyhdf.SD import SD, SDC
from scipy.interpolate import interp1d
import numpy as np
import matplotlib.ticker as plticker

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
    iscan = 600
    icross = 1000

    inRoot = '../outputs'
    lblFile = '{}/hyperTest_LBL_Thuillier_g5nr.nc4'.format(inRoot)
    inRoot = '/nobackup/PACE/LevelC2/Y2006/M03/D24/v2.0'
    ckFile = '{}/OCI2020084005000.L1B_PACE.morel_noaerosol_nosleave.nc'.format(inRoot)
    inRoot = '/nobackup/PACE/LevelC/Y2006/M03/D24/v2.0'
    ckRefFile = '{}/pace-g5nr.lc.vlidort.rayleigh_sleavefree.20060324_005000.nc4'.format(inRoot)
    rsrf0File = '../F0_rsr_weighted_V0.nc'
    rsrFile = 'OCI_RSR_v0.hdf'
    outFile_comparison = 'rsr_v0/l1b_comparison.pdf'
    outFile_difference = 'rsr_v0/l1b_difference.pdf'
    outFile_RSR        = 'rsr_v0/l1b_RSR.pdf'

    lbl = Dataset(lblFile)
    ck  = Dataset(ckFile)
    ckRef = Dataset(ckRefFile)

    # get LBL reflectance
    toa_lbl = lbl.variables['R'][:]
    wav_lbl = np.array(lbl.variables['channels'][:])

    # get sza
    sza = float(lbl.variables['sza'][:])
    csza = np.cos(np.radians(sza))

    # lbl solar irradiance
    F0 = lbl.variables['solar_irradiance'][:]
    # Convert to W/m2/nm
    F0 = F0*1e-2

    # flip this so goes from low to high wavelength
    toa_lbl = toa_lbl[::-1]
    wav_lbl = wav_lbl[::-1]
    F0      = F0[::-1]

    # get CK radiance
    # already integrated over RSR
    wav_blue = ck.groups['sensor_band_parameters'].variables['blue_wavelength'][:]
    wav_red = ck.groups['sensor_band_parameters'].variables['red_wavelength'][:]
    wav_swir = ck.groups['sensor_band_parameters'].variables['SWIR_wavelength'][:]

    sza = ck.groups['geolocation_data'].variables['solar_zenith'][iscan,icross]
    csza = np.cos(np.radians(sza))

    t_blue = ck.groups['observation_data'].variables['Lt_blue'][:,iscan,icross]
    t_red  = ck.groups['observation_data'].variables['Lt_red'][:,iscan,icross]
    t_swir = ck.groups['observation_data'].variables['Lt_SWIR'][:,iscan,icross]

    # get CK reflectance
#    t_blue = ckRef.variables['ref_blue'][iscan,icross,:]
#    t_red  = ckRef.variables['ref_red'][iscan,icross,:]
#    t_swir = ckRef.variables['ref_SWIR'][iscan,icross,:]


    ck_smooth = np.ma.append(t_blue,t_red)
    ck_smooth = np.ma.append(ck_smooth,t_swir)
    wav_ck = np.append(wav_blue,wav_red)
    wav_ck = np.append(wav_ck,wav_swir)

    # integrate over RSR
    # Read in OCI RSR
    inFile = '../{}'.format(rsrFile)
    rsr, wav_rsr, wav_oci = get_rsr(inFile)
    noci = len(wav_oci)
    rsr_f = interp1d(wav_rsr,rsr,kind='linear',fill_value="extrapolate")

    # smooth lbl
    istart_lbl = 0
    rsr_int = rsr_f(wav_lbl)
    lbl_smooth = np.ma.masked_all(noci)
    for ich in range(istart_lbl,noci):
        norm = np.trapz(rsr_int[ich,:],wav_lbl)
        if norm != 0:
            lbl_smooth[ich] = np.trapz(toa_lbl*rsr_int[ich,:],wav_lbl)/norm

    # smooth F0
    F0_smooth = np.ma.masked_all(noci)
    for ich in range(istart_lbl,noci):
        norm = np.trapz(rsr_int[ich,:],wav_lbl)
        if norm != 0:
            F0_smooth[ich] = np.trapz(F0*rsr_int[ich,:],wav_lbl)/norm

    # put in same order as L1B
    # realias 1038 to 1040
    i = wav_oci == 1038
    wav_oci[i] = 1040
    lbl_l1b = np.ma.zeros(ck_smooth.shape)
    F0_l1b  = np.ma.zeros(ck_smooth.shape)
    for ich, ch in enumerate(wav_oci):
        for i,lch in enumerate(wav_ck):
            if float(ch) == lch:
                lbl_l1b[i] = lbl_smooth[ich]
                F0_l1b[i]  = F0_smooth[ich]
    lbl_l1b.mask = ck_smooth.mask
    F0_l1b.mask  = ck_smooth.mask

    # convert radiance to reflectance
    nc = Dataset(rsrf0File)
    F0_rsr = nc.variables['F0_rsr'][:]*1e-2
    F0 = np.ma.zeros(ck_smooth.shape)
    for ich, ch in enumerate(wav_oci):
        for i,lch in enumerate(wav_ck):
            if float(ch) == lch:
                F0[i] = F0_rsr[ich]    
    ck_ref = np.pi*ck_smooth/(csza*F0)
    lbl_ref = np.pi*lbl_l1b/(csza*F0_l1b)
    sys.exit() 
    # ------
    # plotting part comparison
    # -------
    loc = plticker.MultipleLocator(base=50.0)
    xlim = [490,wav_lbl.max()]
    # create figure
    fig = plt.figure()

    ax = fig.add_subplot(2,1,1)

    # lbl
    ax.plot(wav_lbl,toa_lbl,label='LBL')
    # ck
    ax.plot(wav_ck,toa_ck,'o',label='C-K',markersize=1)

    # formatting
    ax.legend()
    ax.set_xlabel('wavelength [nm]')
    ax.set_ylabel('TOA Reflectance')
    ax.set_title('Original')
    ax.set_xlim(xlim)


    # ---- smoothed ----
    ax = fig.add_subplot(2,1,2)

    # lbl
    ax.plot(wav_oci,lbl_smooth,'o',label='LBL',markersize=2)
    # ck
    ax.plot(wav_oci,ck_smooth,'o',label='C-K',markersize=1)

    # formatting
    ax.set_xlabel('wavelength [nm]')
    ax.set_ylabel('TOA Reflectance')
    ax.set_title('Smooth')
    ax.set_xlim(xlim)

    plt.tight_layout()
    plt.savefig(outFile_comparison,bbox_inches='tight')
#    plt.show()
    plt.close()

    # ------
    # plotting part RSR
    # -------

    # create figure
    fig = plt.figure()
    rsr_int = rsr_f(wav_lbl)
    ax = fig.add_subplot(2,1,1)
    for ich in range(noci):
        ax.plot(wav_lbl,rsr_int[ich,:],'k-')

    #formatting
    ax.set_xlabel('wavelength [nm]')
    ax.set_ylabel('OCI Spectral Response Function')
    ax.set_xlim((550,900))
    plt.tight_layout()
    plt.savefig(outFile_RSR,bbox_inches='tight')
#    plt.show()
    plt.close()

    # -------
    # plotting part difference
    # --------

    # create figure
    fig = plt.figure()

    # smooth
    ax = fig.add_subplot(2,1,1)
    diff = 100.*(ck_smooth-lbl_smooth)/lbl_smooth
    ax.plot(wav_oci,diff,'o',label='smooth',markersize=1)
    ax.plot(xlim,[-0.5,-0.5],'k:',linewidth=0.5)
    ax.plot(xlim,[0.5,0.5],'k:',linewidth=0.5)

    # formatting
    ax.set_xlabel('wavelength [nm]')
    ax.set_ylabel('CK - LBL [%]')
    ax.set_title('Smooth')
    ax.set_xlim(xlim)
    ax.yaxis.grid()
    ax.xaxis.set_minor_locator(loc)
    plt.tight_layout()
    plt.savefig(outFile_difference,bbox_inches='tight')
#    plt.show()
    plt.close()
