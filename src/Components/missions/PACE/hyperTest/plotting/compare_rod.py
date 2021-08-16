#!/usr/bin/env python
"""
read and plot o3 spectra from text file
"""

import os
import sys
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import numpy as np
import matplotlib.ticker as plticker

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


#------------------------------------ M A I N ------------------------------------

if __name__ == "__main__":
    outfile_xsec    = 'rod_plots/rod_compare.pdf'
    outFile_difference = 'rod_plots/rod_difference.pdf'

    inRoot = "../"
    inFile = '{}/rod_rsr_weighted_V0.nc'.format(inRoot)

    nc = Dataset(inFile)
    wav_lbl = nc.variables['wav_lbl'][:]
    rod_lbl = nc.variables['rod'][:]
    wav_rsr = nc.variables['wav_oci'][:]
    rod_rsr = nc.variables['rod_rsr'][:]
    nc.close()

    #read OCI table
    wav,rod,depol = read_ROD_table()

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
