#! /usr/bin/env python3
import os, sys
import matplotlib.pyplot as plt
import numpy as np
from glob import glob
import xarray as xr
import matplotlib.pyplot as plt

if __name__ == '__main__':
    row = 85
    col = 220
    rootDir = '/discover/nobackup/projects/gmao/osse2/pub/c1440_NR/OBS/SBG/SSD650'
    vFiles = sorted(glob(rootDir + '/LevelC/Y2006/M01/D16/ssd650-sbg-g5nr.lc.vlidort.20060116_173[0-4]z.nc4'))
    brdfFiles = sorted(glob(rootDir + '/LevelB/Y2006/M01/D16/ssd650-g5nr.lb2.ames_brdf.20060116_173[0-4]z.nc4'))
    aerFiles = sorted(glob(rootDir + '/LevelB/Y2006/M01/D16/ssd650-g5nr.lb2.aer_Nv.20060116_173[0-4]z.nc4'))

    vlid = xr.open_mfdataset(vFiles,chunks='auto')
    brdf = xr.open_mfdataset(brdfFiles,chunks='auto')
    aerm = xr.open_mfdataset(aerFiles,chunks='auto')

    channels = brdf.nwav.values*1e3


    # plot the surface reflectance
    

    riso     = brdf.groups['gridded'].variables['Ki'][:,row,col]
    rgeo     = brdf.groups['gridded'].variables['Kg'][:,row,col]
    rvol     = brdf.groups['gridded'].variables['Kv'][:,row,col]


    vlidort = Dataset(vlidortFile)
    toa = vlidort.variables['toa_reflectance'][:]
    sfc = vlidort.variables['surf_reflectance'][:]

    vlidort.close()

    fig = plt.figure(figsize=(10,8)) 
    ax = fig.add_subplot(1,1,1)

    ax.plot(channels,riso,label='RTLS Riso')
    ax.plot(channels,rgeo,label='RTLS Rgeo')
    ax.plot(channels,rvol,label='RTLS Rvol')

    ax.plot(channels,toa,label='TOA Reflecntace')
    ax.plot(channels,sfc,label='Surfance Reflectance')

    ax.set_title('SBG Simulation 12:30 LT Overpass Nadir View')
    ax.set_ylabel('RTLS Parameters or Reflectance [unitless]')
    ax.set_xlabel('Wavelength [nm]')

    plt.legend()
    plt.savefig('sbg-g5nr.h10v04.20060501_20z.png')
#    plt.show()

