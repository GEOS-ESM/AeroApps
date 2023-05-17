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
    wav_oci = hdf.select('centerwavelength')[:]
    hdf.end()

    return rsr.T, wav_rsr, wav_oci
#------------------------------------ M A I N ------------------------------------

if __name__ == "__main__":

    inRoot = '../outputs'
    lblFile = '{}/hyperTest_LBL_Thuillier.nc4'.format(inRoot)
    ckFile  = '{}/hyperTest_CK_Thuillier.nc4'.format(inRoot)
    rsrFile = 'OCI_RSR_2.5nm_v5.hdf'
    outFile_comparison = 'rsr_v5/comparison.pdf'
    outFile_difference = 'rsr_v5/difference.pdf'
    outFile_RSR        = 'rsr_v5/RSR.pdf'

    lbl = Dataset(lblFile)
    ck  = Dataset(ckFile)

    toa_lbl = lbl.variables['R'][:]
    wav_lbl = np.array(lbl.variables['channels'][:])
    # reverse these so they are going from small to large wavelength
    wav_lbl = wav_lbl[::-1]
    toa_lbl = toa_lbl[::-1]

    toa_ck = ck.variables['R'][:]
    wav_ck = np.array(ck.variables['channels'][:])

    # append lbl values to CK 
    I = wav_lbl < wav_ck.min()
    toa_ck = np.append(toa_lbl[I],toa_ck)
    wav_ck = np.append(wav_lbl[I],wav_ck)


    # integrate over RSR
    # Read in OCI RSR
    inFile = '../{}'.format(rsrFile)
    rsr, wav_rsr, wav_oci = get_rsr(inFile)
    noci = len(wav_oci)
    rsr_f = interp1d(wav_rsr,rsr,kind='linear',fill_value="extrapolate")

    # smooth lbl
    rsr_int = rsr_f(wav_lbl)
    lbl_smooth = np.ma.masked_all(noci)
    for ich in range(noci):
        norm = np.trapz(rsr_int[ich,:],wav_lbl)
        if norm != 0:
            lbl_smooth[ich] = np.trapz(toa_lbl*rsr_int[ich,:],wav_lbl)/norm

    # smooth ck
    istart_ck = 103
    ck_f = interp1d(wav_ck,toa_ck,kind='linear',fill_value='extrapolate')
    toa_ck_int = ck_f(wav_rsr)
    ck_smooth  = np.ma.masked_all(noci)
    for ich in range(istart_ck,noci):
        norm = np.trapz(rsr[ich,:],wav_rsr)
        if norm != 0:
            ck_smooth[ich] = np.trapz(toa_ck_int*rsr[ich,:],wav_rsr)/norm
   
    # ------
    # plotting part comparison
    # -------
    loc = plticker.MultipleLocator(base=50.0)
    xlim = [490,wav_lbl.max()]
    # create figure
    fig = plt.figure()

    # ---- 500+  nm ----
    ax = fig.add_subplot(2,1,1)

    # lbl
    ax.semilogy(wav_lbl,toa_lbl,label='LBL')
    # ck
    ax.semilogy(wav_ck,toa_ck,'o',label='C-K',markersize=1)

    # formatting
    ax.legend()
    ax.set_xlabel('wavelength [nm]')
    ax.set_ylabel('TOA Reflectance')
    ax.set_title('Original')
    ax.set_xlim(xlim)

    # ---- 500+  nm ----
    ax = fig.add_subplot(2,1,2)

    # lbl    
    ax.semilogy(wav_oci,lbl_smooth,'o',label='LBL',markersize=2)
    # ck
    ax.semilogy(wav_oci,ck_smooth,'o',label='C-K',markersize=1)

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
