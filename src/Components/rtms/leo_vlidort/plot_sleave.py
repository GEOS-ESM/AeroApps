#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from matplotlib import cm, colors, rc
from scipy import interpolate

def readRRS(rrsFile,channels):
    """
    sun normalized water leaving radiance
    divided by extraterrestrial solar irradiance
    this has to be multiplied by cos(SZA) to get
    the VLIDORT SL_ISOTROPIC value
    valid for 400-450 nm
    """

    nc = Dataset(rrsFile)
    rrs = nc.variables['rrs'][0,:,900,600]
    wav = nc.variables['wavelength'][:]
    f = interpolate.interp1d(wav,rrs)
    RRS = f(channels)
    i = RRS<0
    RRS[i] = 0.0

    return RRS



if __name__ == '__main__':
    
    code    = 'graspConfig_12_Osku_DrySU_V1_hyperspectral'
    albedo  = 'OCIGissCX_NOBM_CLOUD'
    channel = '350d00_700d00_10'
    inDir = 'self_benchmark_2p8p2/{}/benchmark_rayleigh+simple_aerosol_{}_sleave_adjust/'.format(code,albedo)
    inFile = inDir + 'calipso-g5nr.vlidort.vector.{}.{}.nc4'.format(albedo,channel)
    nc = Dataset(inFile)

    channels   = np.array(nc.variables['channels'][:])
    rrsFile    = '/nobackup/PACE/LevelB/surface/SLEAVE/NOBM/Y2006/M03/D24/pace-g5nr.lb.sleave.20060324_005000.nc4'
    RRS = readRRS(rrsFile,channels)

    RRSa = nc.variables['adjusted_sleave'][:,0,0]
    


    fig = plt.figure()
    ax = plt.subplot(2,1,1)
    ax.plot(channels,RRS,label='input RRS')
    ax.plot(channels,RRSa,label='adjusted RRS')
    ax.set_ylabel('RRS')
    plt.legend()
    ax = plt.subplot(2,1,2)
    ax.plot(channels,RRSa/RRS)
    ax.set_ylabel('Ratio')
    plt.show()



