#!/usr/bin/env python


import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from datetime import datetime, timedelta
from mpl_toolkits.basemap import Basemap
from   dateutil.parser import parse         as isoparser
from glob import glob

if __name__ == "__main__":

    sat = 'SS450'
    inst = 'polar07'
    ialong = 4
    inDir = '/discover/nobackup/projects/gmao/osse2/pub/c1440_NR/OBS/A-CCP/{}/LevelB/Y2006/M06/D20'.format(sat)

    filelist = np.sort(glob(inDir+'/*{}*.nc4'.format(inst)))

    lon, lat, sza, scat, scatss, iso = [],[],[],[],[],[]
    for f in filelist:
        nc = Dataset(f)
        var = nc.variables['scatAngle_ss'][:]
        scatss.append(var)
        var = nc.variables['longitude'][:,ialong,:]
        lon.append(var)
        var = nc.variables['latitude'][:,ialong,:]
        lat.append(var)
        var = nc.variables['scatAngle'][:]
        scat.append(np.max(var,axis=1))
        var = nc.variables['sza'][:,ialong,:]
        sza.append(var)
        var = nc.variables['isotime'][:]
        for v in var:
           iso.append(isoparser(''.join(v)))
        nc.close() 

    lon = np.concatenate(lon,axis=0)
    lat = np.concatenate(lat,axis=0)
    scat = np.concatenate(scat,axis=0)
    scatss = np.concatenate(scatss,axis=0)
    sza  = np.concatenate(sza,axis=0)
    iso = np.array(iso)

    ncross = lon.shape[1]
    mid    = int(ncross*0.5)
    eq = np.abs(lat[:,mid]) < 1
   
    hr = []
    for ll in lon[:,mid]:
        tm = 60.*ll/15.
        hr.append(timedelta(minutes=int(tm)))
    hr = np.array(hr)
    styme = iso + hr

    m = Basemap()
    xx,yy = m(lon,lat)

    scat.mask = sza > 90
    aa = m.contourf(xx,yy,scat,levels=np.linspace(80,180,11))
    plt.colorbar(aa,orientation='horizontal')
    m.drawcoastlines()
    m.drawparallels(range(-80,90,10),labels=np.ones(17))
    plt.title('{} SZA<90'.format(sat))
    plt.savefig('{}_scat90.pdf'.format(sat))
    plt.show()
