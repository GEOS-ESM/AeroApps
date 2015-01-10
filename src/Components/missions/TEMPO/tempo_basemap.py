#!/usr/bin/env python
"""
Reverse engineering of basemap GEO coordinates from lat/lon arrays.
"""
import os

from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap
from matplotlib.pyplot import figure, show, title, savefig
from numpy import arange

#---
def getCoords(m,lon,lat,name):

    # No undefs
    # ---------
    I = (lon<1e14)&(lat<1e14)
    X, Y = m(lon[I],lat[I])

    J = (X<1e30)&(Y<1e30)
    X, Y = X[J], Y[J]

    XA, YA, XB, YB = m.llcrnrx, m.llcrnry, m.urcrnrx, m.urcrnry
    DX, DY = XB-XA, YB-YA
    XC, YC = (XA+XB)/2, (YA+YB)/2

    Xa,Ya,Xb,Yb = (X.min(),Y.min(),X.max(),Y.max())
    xa,ya,xb,yb = (X.min()-XC,Y.min()-YC,X.max()-XC,Y.max()-YC)

    xa_ = (xa)/DX
    xb_ = (xb)/DX
    ya_ = (ya)/DY
    yb_ = (yb)/DY

    BBOX = (Xa, Ya, Xb, Yb)
    Bbox = (xa,ya,xb,yb)
    bbox = (xa_,ya_,xb_,yb_)

    print
    print name
    print 'Native     Bounding box: ', BBOX 
    print 'Recentered Bounding box: ', Bbox 
    print 'Normalized Bounding box: ', bbox 

    return (X,Y,BBOX)

#---
if __name__ == '__main__':

    calculon = '/nobackup/TEMPO/LevelG/invariant/tempo.lg1.invariant.nc4'
    nccs = '/home/dasilva/silo/TEMPO/LevelG/invariant/tempo.lg1.invariant.nc4'

    # Get coordinates
    # ---------------
    if os.path.exists(calculon): 
        nc = Dataset(calculon)
    else:
        nc = Dataset(nccs)

    # Read and Tighten the E-W domain
    # -------------------------------
    clon = nc.variables[u'clon'][:,124:1374]
    clat = nc.variables[u'clat'][:,124:1374]
    elon = nc.variables[u'elon'][:,124:1374]
    elat = nc.variables[u'elat'][:,124:1374]

    mx, my = elon.shape
    nx, ny = clon.shape

    # GEO disk parked at 100W
    # -----------------------
    lon_0 = -100.
    m =  Basemap(projection='geos',lon_0=lon_0,resolution=None,
                 rsphere=(6378137.00,6356752.3142),
                 satellite_height = 35785831.0)

    cX, cY, cBOX = getCoords(m,clon,clat,'CENTER Coordinates')
    eX, eY, eBOX = getCoords(m,elon,elat,'EDGE   Coordinates')


    # Plot Map
    # --------
    m =  Basemap(projection='geos',lon_0=lon_0,resolution='l',
                 llcrnrx=eBOX[0],
                 llcrnry=eBOX[1],
                 urcrnrx=eBOX[2],
                 urcrnry=eBOX[3],
                 rsphere=(6378137.00,6356752.3142),
                 satellite_height = 35785831.0)
    

    figure()
    m.drawcoastlines()
    m.fillcontinents(color='coral',lake_color='aqua')
    m.drawcountries()
    m.drawstates()
    m.drawparallels(arange(-90.,120.,5.),labels=[0,0,0,0])
    m.drawmeridians(arange(0.,360.,5.),labels=[0,0,0,0])
    m.drawmapboundary(fill_color='aqua')

    title('TEMPO Geostationary Sector')
    savefig('tempo_basemap.png')
    #show()
