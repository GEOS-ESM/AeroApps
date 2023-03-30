#!/usr/bin/env python3


import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from datetime import datetime, timedelta
from mpl_toolkits.basemap import Basemap
from   dateutil.parser import parse         as isoparser
from glob import glob

if __name__ == "__main__":

    sat = 'AQUA'
    inst = 'polder'
    ialong = 74
    inDir = '/discover/nobackup/projects/gmao/osse2/pub/c1440_NR/OBS/A-CCP/{}/LevelB/Y2006/M06/D20'.format(sat)

    filelist = np.sort(glob(inDir+'/*{}*.nc4'.format(inst)))

    lon, lat, sza, saa , scat, iso = [],[],[],[],[],[]
    for f in filelist:
        nc = Dataset(f)
        var = nc.variables['longitude'][:,ialong,:]
        lon.append(var)
        var = nc.variables['latitude'][:,ialong,:]
        lat.append(var)
        var = nc.variables['scatAngle'][:,ialong,:]
        scat.append(var)
        var = nc.variables['sza'][:,ialong,:]
        sza.append(var)
        var = nc.variables['saa'][:,ialong,:]
        saa.append(var)
        var = nc.variables['isotime'][:]
        for v in var:
           iso.append(isoparser(''.join(v)))
        nc.close() 

    lon = np.concatenate(lon,axis=0)
    lat = np.concatenate(lat,axis=0)
    scat = np.concatenate(scat,axis=0)
    sza  = np.concatenate(sza,axis=0)
    saa  = np.concatenate(saa,axis=0)
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

    hr = []
    for tyme in iso:
        hh = tyme.hour + tyme.minute/60.
        hr.append(hh)

    hr = np.ma.array(hr)


    m = Basemap()
    xx,yy = m(lon,lat)

    I  = sza > 90
    sza.mask = I
    saa.mask = I
    hr.mask  = I

    doplot = True
    if doplot:
        aa = m.contourf(xx,yy,saa,levels=np.linspace(0,360,37),cmap=plt.get_cmap('jet'))
        plt.colorbar(aa,orientation='horizontal')
        m.drawcoastlines()
        m.drawparallels(list(range(-80,90,10)),labels=np.ones(17))
        m.drawmeridians(list(range(-170,180,10)),labels=np.ones(len(list(range(-170,180,10)))))
        plt.title('{} SAA'.format(sat))
        plt.savefig('{}_saa90nadir.pdf'.format(sat))
        plt.show()

        aa = m.contourf(xx,yy,sza,levels=np.linspace(0,90,19),cmap=plt.get_cmap('jet'))
        plt.colorbar(aa,orientation='horizontal')
        m.drawcoastlines()
        m.drawparallels(list(range(-80,90,10)),labels=np.ones(17))
        m.drawmeridians(list(range(-170,180,10)),labels=np.ones(len(list(range(-170,180,10)))))
        plt.title('{} SZA'.format(sat))
        plt.savefig('{}_sza90nadir.pdf'.format(sat))
        plt.show()


        aa = m.scatter(xx[:,mid],yy[:,mid],c=hr,cmap=plt.get_cmap('jet'))
        plt.colorbar(aa,orientation='horizontal')
        m.drawcoastlines()
        m.drawparallels(list(range(-80,90,10)),labels=np.ones(17))
        m.drawmeridians(list(range(-170,180,10)),labels=np.ones(len(list(range(-170,180,10)))))
        plt.title('{} UTC'.format(sat))
        plt.savefig('{}_utc90nadir.pdf'.format(sat))
        plt.show()

