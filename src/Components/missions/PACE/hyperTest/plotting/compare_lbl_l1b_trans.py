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
    iscan  = 600
    icross = 1000
    inRoot = '../outputs'
    lblFile = '{}/hyperTest_LBL_Thuillier_g5nr.nc4'.format(inRoot)
    inRoot = '/nobackup/PACE/LevelC2/Y2006/M03/D24/v2.0'
    geoFile  = '{}/OCI2020084005000.L1B_PACE.nc'.format(inRoot)
    inRoot = '/nobackup/PACE/LevelC/Y2006/M03/D24/v2.0'
    ckFile = '{}/pace-g5nr.lc.gas.20060324_005000.nc4'.format(inRoot)
    rsrFile = 'OCI_RSR_v0.hdf'
    outFile_comparison = 'rsr_v0/l1b_trans_comparison.pdf'
    outFile_difference = 'rsr_v0/l1b_trans_difference.pdf'
    outFile_RSR        = 'rsr_v0/l1b_trans_RSR.pdf'

    lbl = Dataset(lblFile)
    ck  = Dataset(ckFile)
    geo = Dataset(geoFile)

    # get angle
    vza = geo.groups['geolocation_data'].variables['sensor_zenith'][iscan,icross]
    cvza = np.cos(np.radians(vza))

    # get LBL trans
    alpha_lbl = lbl.variables['ROD'][:]
    toa_lbl = np.exp(-1.*alpha_lbl/cvza)
    wav_lbl = np.array(lbl.variables['channels'][:])

    # lbl solar irradiance
    F0 = lbl.variables['solar_irradiance'][:]

    # flip this so goes from low to high wavelength
    toa_lbl = toa_lbl[::-1]
    wav_lbl = wav_lbl[::-1]
    F0      = F0[::-1]

    # get CK transmittance
    # already integrated over RSR
    wav_blue = ck.variables['blue_wavelength'][:]
    wav_red = ck.variables['red_wavelength'][:]
    wav_swir = ck.variables['SWIR_wavelength'][:]

    t_blue = ck.variables['TRANS_H2O_blue'][iscan,icross,:]
    t_red  = ck.variables['TRANS_H2O_red'][iscan,icross,:]
    t_swir = ck.variables['TRANS_H2O_SWIR'][iscan,icross,:]

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
    istart_lbl = 1
    rsr_int = rsr_f(wav_lbl)
    lbl_smooth = np.ma.masked_all(noci)
    for ich in range(istart_lbl,noci):
        norm = np.trapz(rsr_int[ich,:]*F0,wav_lbl)
        if norm != 0:
            lbl_smooth[ich] = np.trapz(F0*toa_lbl*rsr_int[ich,:],wav_lbl)/norm


    # put in same order as L1B
    # realias 1038 to 1040
    i = wav_oci == 1038
    wav_oci[i] = 1040
    lbl_l1b = np.ma.zeros(ck_smooth.shape)
    for ich, ch in enumerate(wav_oci):
        for i,lch in enumerate(wav_ck):
            if float(ch) == lch:
                lbl_l1b[i] = lbl_smooth[ich]
    lbl_l1b.mask = ck_smooth.mask
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
