#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from matplotlib import cm, colors, rc

if __name__ == '__main__':

    albedo  = 'OCIGissCX_NOBM_CLOUD'
    channel = '400'
    inDir = 'self_benchmark_2p8p1/graspConfig_12_Osku_DrySU_V1/benchmark_rayleigh+simple_aerosol_{}/'.format(albedo)
    inFile = inDir + 'calipso-g5nr.vlidort.vector.{}.{}d00.nc4'.format(albedo,channel)
    nc = Dataset(inFile)
    inDir = 'self_benchmark_2p8p2/graspConfig_12_Osku_DrySU_V1/benchmark_rayleigh+simple_aerosol_{}/'.format(albedo)
    inFile = inDir + 'calipso-g5nr.vlidort.vector.{}.{}d00.nc4'.format(albedo,channel)
    ncn = Dataset(inFile)

    I = nc.variables['I'][:]
    In = ncn.variables['I'][:]
    diff = I - In

    vmax = np.max([I.max(),In.max()])
    vmin = np.min([I.min(),In.min()])

    vza = nc.variables['sensor_zenith'][:]
    vaa = nc.variables['sensor_azimuth'][:]

    r, theta = np.meshgrid(vza,vaa)
    r = r.T
    theta = theta.T    
    theta = np.radians(theta)

    fig = plt.figure()
    ax = plt.subplot(2,1,1,polar=True)
    ax.contourf(theta,r,I,vmin=vmin,vmax=vmax)
    ax = plt.subplot(2,1,2,polar=True)
    ax.contourf(theta,r,In,vmin=vmin,vmax=vmax)
    plt.show()
    ax = plt.subplot(1,1,1,polar=True)
    cax = ax.contourf(theta,r,diff)
    cbar = fig.colorbar(cax,ax=ax,orientation='horizontal')
    plt.show()



