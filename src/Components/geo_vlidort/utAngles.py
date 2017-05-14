#!/usr/bin/env python
"""
   Unit testing for Geostationary angles.
"""
from datetime        import datetime, timedelta

from GeoAngles_ import satangles
from numpy      import ones, arange, array, tile

def scanTimes(nEW,nNS,tBeg):
    """
    Calculate TEMPO scantimes
    """
    tScan = array([tBeg + i * timedelta(seconds=2.85) for i in arange(nEW) ])
    return tile(tScan,(nNS,1))

if __name__ == "__main__":

    yr=2014 * ones(3,dtype='int32')
    mo=8 * ones(3,dtype='int32')
    day=1 * ones(3,dtype='int32')
    hr=18 * ones(3,dtype='int32')
    min=0 * ones(3,dtype='int32')
    sec=0 * ones(3,dtype='int32')
    latp=36.5 * ones(3)
    lonp=-100.0 * ones(3)
    lonss=-100.0

    vaa, vza, saa, sza = satangles(yr,mo,day,hr,min,sec,latp,lonp,lonss)

    print "Angles:"
    print "Viewing azimuth: ", vaa
    print "Viewing  zenith: ", vza
    print "Solar   azimuth: ", saa
    print "Solar    zenith: ", sza

    nEW,nNS = 5, 3

    tBeg = datetime(2014,11,20,11)
    t = scanTimes(nEW, nNS,tBeg)

