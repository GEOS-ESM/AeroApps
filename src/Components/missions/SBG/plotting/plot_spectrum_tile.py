#! /usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset


if __name__ == '__main__':
    row = 85
    col = 220
    vlidortFile = '/discover/nobackup/projects/gmao/osse2/pub/c1440_NR/OBS/SBG/LevelC/Y2006/M05/sbg-g5nr.vlidort.AMES_BRDF.20060501_20z.nc4'
    brdfFile = '/discover/nobackup/projects/gmao/osse2/pub/c1440_NR/OBS/SBG/LevelB/surface/BRDF/ames-brdf/version2.0/2019/h10v04/sbg_mdl_h10v04_2019121_brdf.nc4'

    brdf = Dataset(brdfFile)

    channels = brdf.groups['ancillary'].variables['wavelength'][:]*1e3
    riso     = brdf.groups['gridded'].variables['Ki'][:,row,col]
    rgeo     = brdf.groups['gridded'].variables['Kg'][:,row,col]
    rvol     = brdf.groups['gridded'].variables['Kv'][:,row,col]

    brdf.close()

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

