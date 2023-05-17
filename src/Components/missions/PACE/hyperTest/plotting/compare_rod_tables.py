#!/usr/bin/env python3
"""
read and plot o3 spectra from text file
"""

import os
import sys
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import numpy as np
import matplotlib.ticker as plticker
from pyhdf.SD import SD, SDC

def read_ROD_table():
    f = open('../../../../leo_vlidort/amir/OCI_ROD_Table_adjusted.txt','r')
    f.readline() #header
    wav = f.readline().split()
    f.readline() #wav center
    f.readline() #wav width
    f.readline() #F0
    rod = f.readline().split()
    depol = f.readline().split()

    f.close()

    rod = rod[2:]
    rod = np.array(rod).astype('float')

    depol = depol[2:]
    depol = np.array(depol).astype('float')

    wav = wav[3:]
    wav = np.array(wav).astype('float')

    return wav,rod,depol


def read_highres_table():
    f = open('../../oci_tables/rayleigh_bodhaine.txt')

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
    #read OCI tables
    wav,rod,depol = read_ROD_table()
    hwav,hrod,hdepol = read_highres_table()

    rsrFile = 'OCI_RSR_v0.hdf'
    # Read in OCI RSR
    inFile = '../{}'.format(rsrFile)
    rsr, wav_rsr, wav_oci = get_rsr(inFile)
    noci = len(wav_oci)
    rsr_f = interp1d(wav_rsr,rsr,kind='linear',fill_value="extrapolate")
 
    i = (hwav<=wav_rsr.max()) & (hwav>= wav_rsr.min())
    hwav = hwav[i]
    hrod = hrod[i]
    hdepol = hdepol[i]
    

    # smooth hired
    rsr_int = rsr_f(hwav)
    rod_smooth = np.zeros(noci)
    for ich in range(noci):
        norm = np.trapz(rsr_int[ich,:],hwav)
        rod_smooth[ich] = np.trapz(hrod*rsr_int[ich,:],hwav)/norm

    sys.exit()
    # xsec
    loc = plticker.MultipleLocator(base=50.0)
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.semilogy(wav_lbl,rod_lbl,label='LBL')
    ax.semilogy(wav_rsr,rod_rsr,label='RSR weighted')
    ax.semilogy(wav,rod,label='OCI Table')

    ax.legend()
    ax.set_ylabel('Rayleight Optical Depth')
    ax.set_xlabel('wavelength [nm]')
    ax.xaxis.set_minor_locator(loc)
    plt.savefig(outfile_xsec,bbox_inches='tight')
#    plt.show()
    plt.close()


    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    diff = 100.*(rod_rsr - rod)/rod
    ax.plot(wav,diff)
    ax.set_ylabel('% Difference')
    ax.yaxis.grid()
    ax.xaxis.set_minor_locator(loc)
    plt.tight_layout()
    plt.savefig(outFile_difference,bbox_inches='tight')
#    plt.show()
    plt.close()
