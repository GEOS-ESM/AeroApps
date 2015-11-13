#!/usr/bin/env python
"""
Make the lg1 file for GEMS
"""
import numpy as np
from math import radians
from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
from matplotlib import cm, colors

def set_basemap(projection,clon,clat,elon,elat,lon_0=128.2,satellite_height=35785831.0):
    
  # Default is GEMS
  # GEO disk parked at 128.2E
  # -----------------------
  m =  Basemap(projection=projection,lon_0=lon_0,resolution=None,
       rsphere=(6378137.00,6356752.3142),
       satellite_height = satellite_height)

  cX, cY, cBOX,cbox = getCoords(m,clon,clat,'CENTER Coordinates')
  eX, eY, eBOX,ebox = getCoords(m,elon,elat,'EDGE   Coordinates')

  # Plot Map
  # --------
  m =  Basemap(projection=projection,lon_0=lon_0,resolution='l',                 
       rsphere=(6378137.00,6356752.3142),
       satellite_height = satellite_height,
       llcrnrx=ebox[0],
       llcrnry=ebox[1],
       urcrnrx=ebox[2],
       urcrnry=ebox[3])    

  return m

def map_(fig,m,data,cmap,minval,maxval,title,norm=None,format=None,elon=None,elat=None,
         fill_color=None,bluemarble=None):
    if (m.projection == 'lcc'):
        if (bluemarble is not None):
            m.warpimage(image=bluemarble)
        else:
            m.bluemarble()
        im = m.pcolormesh(elon,elat,data,latlon=True,cmap=cmap,vmin=minval,vmax=maxval,norm=norm)
        im.set_rasterized(True)
    elif (m.projection == 'geos'):
        im = m.imshow(data,cmap=cmap,vmin=minval,vmax=maxval,norm=norm)
    m.drawcoastlines(color='white')
    m.drawcountries(color='white')
    m.drawstates(color='white')

    if (m.projection == 'geos'):
        if (fill_color is None):
            fill_color = 'k' #'0.90'
        m.drawmapboundary(color='k',fill_color=fill_color,linewidth=2.0)
    
    plt.title(title, fontsize=20)

    return im

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

    return (X,Y,BBOX,Bbox)
#---
if __name__ == '__main__':

	outfile    = 'gems.lg1.invariant.nc4'
	projection = 'geos'
	lon_0      = 128.2
	satellite_height = 35785831.0

	centerlon = 114.1  #degrees E
	centerlat = 17.04  #degrees N
	fov       = 70  #degree of view East-West
	latmax = 45
	latmin = -5
	nEW   = 1250/2   #half the TEMPO resoultion
	nNS   = 2000

  #-------
  ##  End of User Input
  #-------
	m =  Basemap(projection=projection,lon_0=lon_0,resolution=None,
       rsphere=(6378137.00,6356752.3142),
       satellite_height = satellite_height)


	# convert FOV to meters on surface
	rearth = np.mean([6378137.00,6356752.3142])
	DX   = rearth * radians(fov)

	# Get sector center in map projection coordinates
	cX, cY = m(centerlon,centerlat)

	maxX = cX + 0.5*DX
	minX = cX - 0.5*DX

	# Get latmax and latmin in map projection coordinates
	cXmin, minY = m(centerlon,latmin)
	cXmax, maxY = m(centerlon,latmax)

	DY = maxY - minY

	dx = DX/nEW
	dy = DY/nNS

	clon = np.ma.masked_all([nNS,nEW])
	clat = np.ma.masked_all([nNS,nEW])
	elon = np.ma.masked_all([nNS+1,nEW+1])
	elat = np.ma.masked_all([nNS+1,nEW+1])	
	for i in np.arange(nNS+1):
		for j in np.arange(nEW+1):
			elon[i,j], elat[i,j] = m(minX + j*dx, minY + i*dy, inverse=True)

	for i in np.arange(nNS):
		for j in np.arange(nEW):			
			clon[i,j], clat[i,j] = m(minX + j*dx + 0.5*dx, minY + i*dy + 0.5*dy, inverse=True)

	Data = np.ones([nNS,nEW])
	Data[np.where(clon == 1e30)] = np.ma.masked
	Data[np.where(clat == 1e30)] = np.ma.masked
		
	cmap = cm.jet
	cmap.set_bad(color='w',alpha=0)

	m = set_basemap(projection,clon,clat,elon,elat,lon_0=128.2,satellite_height=35785831.0)
	fig = plt.figure()
	im = map_(fig,m,Data,cmap,0,2,'',elon=elon,elat=elat)
	
	