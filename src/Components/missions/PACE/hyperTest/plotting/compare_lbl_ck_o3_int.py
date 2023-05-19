#!/usr/bin/env python3

"""
Compare LBL and C-K simulations
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

    inRoot = '../outputs'
    lblFile = '{}/hyperTest_LBL_Thuillier_o3.nc4'.format(inRoot)
    ckFile  = '{}/hyperTest_CK_Thuillier_o3.nc4'.format(inRoot)
    uvFile  = '{}/hyperTest_UV_Thuillier_o3.nc4'.format(inRoot)
    rsrFile = 'OCI_RSR_v0.hdf'
    outFile_comparison = 'o3_rsr_v0/comparison.pdf'
    outFile_difference = 'o3_rsr_v0/difference.pdf'
    outFile_RSR        = 'o3_rsr_v0/RSR.pdf'

    lbl = Dataset(lblFile)
    ck  = Dataset(ckFile)
    uv  = Dataset(uvFile)

    toa_lbl = lbl.variables['R'][:]
    wav_lbl = np.array(lbl.variables['channels'][:])

    # flip this so goes from low to high wavelength
    toa_lbl = toa_lbl[::-1]
    wav_lbl = wav_lbl[::-1]

    toa_ck = ck.variables['R'][:]
    wav_ck = np.array(ck.variables['channels'][:])

    toa_uv = uv.variables['R'][:]
    wav_uv = np.array(uv.variables['channels'][:])

    # regrid CK so a consistent resolution
    wav_ck_grid = np.arange(wav_ck.min(),wav_ck.max()+0.5,0.5)
    toa_ck_f = interp1d(wav_ck,toa_ck,kind='linear')
    toa_ck_grid = toa_ck_f(wav_ck_grid)
    toa_ck = toa_ck_grid
    wav_ck = wav_ck_grid
    
    # integrate over RSR
    # Read in OCI RSR
    inFile = '../{}'.format(rsrFile)
    rsr, wav_rsr, wav_oci = get_rsr(inFile)
    noci = len(wav_oci)
    rsr_f = interp1d(wav_rsr,rsr,kind='linear',fill_value=0.0,bounds_error=False)
    
    # smooth lbl
    rsr_int = rsr_f(wav_lbl)
    lbl_smooth = np.ma.masked_all(noci)
    for ich in range(noci):
        norm = np.trapz(rsr_int[ich,:],wav_lbl)
        if norm != 0:
            lbl_smooth[ich] = np.trapz(toa_lbl*rsr_int[ich,:],wav_lbl)/norm

    # smooth ck
    # but only for OCI channels covered by CK bins
    #
    # ck min 575, i = 53
    # for 2.5 ck min is 570 i = 103
    istart = 53
    rsr_int = rsr_f(wav_ck)
    ck_smooth  = np.ma.masked_all(noci-istart)
    for ich in range(istart,noci):
        norm = np.trapz(rsr_int[ich,:],wav_ck)
        if norm != 0:
            ck_smooth[ich-istart] = np.trapz(toa_ck*rsr_int[ich,:],wav_ck)/norm


    # smooth UV wavelengths less than 330
    ismooth = np.argmax(wav_uv ==  342.75) + 1
    nint = np.sum(ismooth)
    nint = ismooth
    nsmooth = 4

    # interpolate to rsr wavelengths
    #uv_f = interp1d(wav_uv[:ismooth],toa_uv[:ismooth],kind='linear',fill_value='extrapolate',bounds_error=False)
    #uv_int = uv_f(wav_rsr)
    rsr_int = rsr_f(wav_uv[:ismooth])
    uv_smooth = np.ma.masked_all(nsmooth)
    wav_uv_smooth = wav_uv[:ismooth]
    for ich in range(nsmooth):
        norm = np.trapz(rsr_int[ich,:],wav_uv_smooth)
        if norm != 0:
            uv_smooth[ich] = np.trapz(toa_uv[:ismooth]*rsr_int[ich,:],wav_uv_smooth)/norm

    uv_smooth = np.append(uv_smooth,toa_uv[nint:])

    # append RSR weighted UV and CK
    ck_smooth = np.append(uv_smooth[:istart],ck_smooth)

    # ------
    # plotting part comparison
    # -------
    xlim = [wav_lbl.min(),wav_lbl.max()]

    # create figure
    fig = plt.figure()

    # ---- 500+  nm ----
    ax = fig.add_subplot(2,1,1)

    # lbl
    ax.semilogy(wav_lbl,toa_lbl,label='LBL')
    # ck
    ax.semilogy(wav_ck,toa_ck,'o',label='C-K',markersize=1)
    # uv
    ax.semilogy(wav_uv[:istart],toa_uv[:istart],label='UV - O3 only',markersize=1)

    # formatting
    ax.legend()
    ax.set_xlabel('wavelength [nm]')
    ax.set_ylabel('TOA Reflectance')
    ax.set_title('Original')
    ax.set_xlim(xlim)


    # ---- Smoothed ----
    ax = fig.add_subplot(2,1,2)

    # lbl
    I = wav_oci >= wav_lbl.min()
    ax.semilogy(wav_oci[I],lbl_smooth[I],'o',label='LBL',markersize=2)
    # ck
    ax.semilogy(wav_oci[I],ck_smooth[I],'o',label='C-K + UV',markersize=1)

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
    ax.set_xlim((300,900))

    plt.tight_layout()
    plt.savefig(outFile_RSR,bbox_inches='tight')
#    plt.show()
    plt.close()

    # -------
    # plotting part difference
    # --------
    loc = plticker.MultipleLocator(base=50.0)
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
    ax.xaxis.set_minor_locator(loc)
    ax.yaxis.grid()

    plt.tight_layout()
    plt.savefig(outFile_difference,bbox_inches='tight')
#    plt.show()
    plt.close()
