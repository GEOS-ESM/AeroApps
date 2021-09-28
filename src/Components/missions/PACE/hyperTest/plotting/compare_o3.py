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

def get_o3_profile(inFile):
    """
    read o3 profile (molecules/cm3) from mcclams.dat
    p = hPa
    t = K
    airdens = molecules/cm3
    """
    f = open(inFile)
    # header
    f.readline()
    f.readline()

    z = []
    p = []
    t = []
    airdens = []
    o3_conc = []
    for l in f:
        aa = np.array(l.split()).astype(float)
        z.append(aa[0])
        p.append(aa[1])
        t.append(aa[2])
        airdens.append(aa[3])
        o3_conc.append(aa[4])

    z = np.array(z)
    p = np.array(p)
    t = np.array(t)
    airdens = np.array(airdens)
    o3_conc = np.array(o3_conc)

    # o3 mixing ratio, unitless
    o3_vmr = o3_conc/airdens

    # get altitude in m
    z = z*1e3


    f.close()
    return z,o3_vmr



#------------------------------------ M A I N ------------------------------------

if __name__ == "__main__":
    outfile_profile = 'o3_plots/profile.pdf'
    outfile_xsec    = 'o3_plots/xsec_compare.pdf'

    inRoot = "../o3_bremen"
    inFile = '{}/xsec_o3_lbl.nc'.format(inRoot)

    nc = Dataset(inFile)
    wav_o3 = nc.variables['wav_o3'][:]
    xsec   = nc.variables['xsec_o3'][:]
    wav_lbl = nc.variables['wav_lbl'][:]
    xsec_lbl = nc.variables['xsec_o3_lbl'][:]
    nc.close()

    inFile = '{}/xsec_o3_ck_bins_V0.nc'.format(inRoot)
    nc = Dataset(inFile)
    wav_ck = nc.variables['wav_ck'][:]
    xsec_ck = nc.variables['xsec_o3_ck'][:]
    nc.close()

    inFile ='{}/xsec_o3_rsr_weighted_V0.nc'.format(inRoot)
    nc = Dataset(inFile)
    wav_rsr = nc.variables['wav_oci'][:]
    xsec_rsr = nc.variables['xsec_o3_rsr'][:]
    nc.close()
    wav_rsr = wav_rsr[:53]
    xsec_rsr = xsec_rsr[:53,:]


    # Read in mcclams ozone mixing ratio profile [unitless]
    inFile = '../mcclams.dat'
    z_o3, o3_vmr = get_o3_profile(inFile)


    # profile
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.plot(o3_vmr*1e6,z_o3)
    ax.set_xlabel('o3 [ppm]')
    ax.set_ylabel('altitude [m]')
    plt.savefig(outfile_profile,bbox_inches='tight')
    plt.close()

    # xsec
    loc = plticker.MultipleLocator(base=50.0)
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.semilogy(wav_o3,xsec[:,0],label='Bremen')
    ax.semilogy(wav_lbl,xsec_lbl[:,0],label='LBL')
    ax.semilogy(wav_ck,xsec_ck[:,0],label='CK')
    ax.semilogy(wav_rsr,xsec_rsr[:,0],label='RSR weighted')

    ax.legend()
    ax.set_ylabel('O3 Xsection [cm2/molecule]')
    ax.set_xlabel('wavelength [nm]')
    ax.xaxis.set_minor_locator(loc)
    plt.savefig(outfile_xsec,bbox_inches='tight')
    plt.close()
