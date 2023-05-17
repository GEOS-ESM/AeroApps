#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from matplotlib import cm, colors, rc

if __name__ == '__main__':

    inFile = '/home/pcastell/workspace/vlidort_benchmark2p8p1/B3_benchmark.nstreams=20.nc4'
    nc = Dataset(inFile)
    inFile = 'self_benchmark_2p8p1/b3/B3_benchmark.nc4'
    ncn = Dataset(inFile)

    I = nc.variables['I'][0,:,0,:,0].T
    In = ncn.variables['I'][:]
    diff = I - In

    vmax = np.max([I.max(),In.max()])
    vmin = np.min([I.min(),In.min()])

    vza = ncn.variables['sensor_zenith'][:]
    vaa = ncn.variables['sensor_azimuth'][:]

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



