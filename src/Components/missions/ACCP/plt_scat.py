#!/usr/bin/env python3


import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from datetime import datetime, timedelta
from mpl_toolkits.basemap import Basemap
from   dateutil.parser import parse         as isoparser
from glob import glob

def scatAngle(raa,sza,vza):
    raar = np.radians(raa)
    szar = np.radians(sza)
    vzar = np.radians(vza)
    cscat = np.cos(szar)*np.cos(vzar) + np.sin(szar)*np.sin(vzar)*np.cos(raar)
    scat = 180. - np.degrees(np.arccos(cscat))

    return scat

if __name__ == "__main__":
    #June 20

    # North  Lat=35 Lon =4
    saaN = 223.
    szaN = 15.

    # South  Lat=-25.3  Lon=17.5 
    saaS = 331.85
    szaS = 54.26

    # sat
    vaaf = 168.
    vaab = 348.

    # vza
    vza = np.arange(90)
    vzap = np.array([6.74,20.26,33.94,48.02,63.15])

    # RAA
    raaNf = saaN - vaaf
    raaNb = vaab - saaN 

    raaSf = saaS - vaaf
    raaSb = vaab - saaS 

    scatNf = scatAngle(raaNf,szaN,vza)
    scatNb = scatAngle(raaNb,szaN,vza)
    scatSf = scatAngle(raaSf,szaS,vza)
    scatSb = scatAngle(raaSb,szaS,vza)
    scatSpf = scatAngle(raaSf,szaS,vzap)
    scatSpb = scatAngle(raaSb,szaS,vzap)

#    plt.plot(vza,scatNf,label='forward')
#    plt.plot(-1.*vza,scatNb,label='backward')
#    plt.plot([szaN,szaN],[0,180],label='SZA')
#    plt.legend()
#    plt.show()
    plt.plot(vza,scatSf,label='forward')
    plt.plot(-1.*vza,scatSb,label='backward')
    plt.plot([szaS,szaS],[0,180],label='SZA')
    plt.plot([-1.*szaS,-1.*szaS],[0,180],'g-')
    plt.plot(vzap,scatSpf,'ko',label='polarimeter07')
    plt.plot(-1.*vzap,scatSpb,'ko')
    plt.plot([-90,90],[165,165],label='165')
    plt.legend()
    plt.savefig('scat_vza.pdf')
    plt.show()
