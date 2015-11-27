#!/usr/bin/env python
"""
Make the lg1 file for GEMS
"""
import numpy as np
from math import radians, degrees
import math
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
	fovEW     = 10.91  #degree of view East-West
	latmax = 45
	latmin = -5
	nEW   = 1250/2   #half the TEMPO resoultion
	nNS   = 1100

  #-------
  ##  End of User Input
  #-------
	m =  Basemap(projection=projection,lon_0=lon_0,resolution=None,
       rsphere=(6378137.00,6356752.3142),
       satellite_height = satellite_height)

	# calculate total geostationary field of view half angle
	rearth = np.mean([6378137.00,6356752.3142])
	Z = rearth + satellite_height

	alpha = math.asin(rearth/Z)
	fovEW = radians(fovEW)

	# EW viewing angles to consider
	# positive is west of nadir, negative east of nadir
	angles = alpha - (fovEW/nEW)*np.arange(nEW+1)

	# convert FOVEW to distance in meters on surface between nadir and the point
	#DX = np.ma.masked_all(nEW+1)
	DX = np.zeros(nEW+1)
	for i, alphaprime in enumerate(abs(angles)):
		alphaprime = abs(alphaprime)
		beta = 0.5*math.pi - alphaprime
		D = Z*math.tan(alphaprime)
		d = math.asin(D*math.sin(beta)/rearth)
		thetaprime = d + beta - 0.5*math.pi
		DX[i]   = rearth * thetaprime
		if angles[i] < 0:
			DX[i] = -1.*DX[i]

	# EW size of each pixel
	dX = abs(DX[0:-1] - DX[1:])
  
  # Get center X position relative to scan sector center 
	#X = np.ma.masked_all(nEW)
	X = np.zeros(nEW)
	for i, deltax in enumerate(dX[0:int(0.5*(nEW+1))]):
		X[i] = -1.*np.sum(dX[i:int(0.5*(nEW+1))]) + 0.5*dX[i]

	for i, deltax in enumerate(dX[int(0.5*(nEW+1)):]):
		X[i+int(0.5*(nEW+1))] = np.sum(dX[int(0.5*(nEW+1)):int(0.5*(nEW+1))+i]) + 0.5*dX[i+int(0.5*(nEW+1))]


	# Calculate N-S field of view alpha (5S to 45N)
	alphaNS_north = math.atan(rearth*math.sin(radians(45))/(satellite_height - rearth*math.cos(radians(45))))
	alphaNS_south = math.atan(rearth*math.sin(radians(5))/(satellite_height - rearth*math.cos(radians(5))))

	fovNS = alphaNS_north + alphaNS_south

	# NS viewing angles to consider
	# positive is north of nadir, negative south of nadir
	angles = alphaNS_north - (fovNS/nNS)*np.arange(nNS+1)

	# convert FOVNS to distance in meters on surface between nadir and the point
	DY = np.zeros(nNS+1)
	for i, alphaprime in enumerate(abs(angles)):
		beta = 0.5*math.pi - alphaprime
		D = Z*math.tan(alphaprime)
		d = math.asin(D*math.sin(beta)/rearth)
		thetaprime = d + beta - 0.5*math.pi
		DY[i]   = rearth * thetaprime
		if angles[i] < 0:
			DY[i] = -1.*DY[i]

	# NS size of each pixel
	dY = abs(DY[0:-1] - DY[1:])
	dY = dY[::-1]     #reverse so it's south to north

  # Get center Y position relative to scan sector center 
	Y = np.zeros(nNS)
	for i, deltay in enumerate(dY[0:int(0.5*(nNS+1))]):
		Y[i] = -1.*np.sum(dY[i:int(0.5*(nNS+1))]) + 0.5*dY[i] 

	for i, deltay in enumerate(dY[int(0.5*(nNS+1)):]):
		Y[i+int(0.5*(nNS+1))] = np.sum(dY[int(0.5*(nNS+1)):int(0.5*(nNS+1))+i]) + 0.5*dY[i+int(0.5*(nNS+1))] 

	# Because nNS is even
	Y = Y - 0.5*dY[int(0.5*(nNS+1))]

	# Get sector center in map projection coordinates
	cX, cY = m(centerlon,centerlat)

	cX = cX + X
	cY = cY + Y
	
	# clon = np.ma.masked_all([nNS,nEW])
	# clat = np.ma.masked_all([nNS,nEW])
	# elon = np.ma.masked_all([nNS+1,nEW+1])
	# elat = np.ma.masked_all([nNS+1,nEW+1])	

	clon = np.zeros([nNS,nEW])
	clat = np.zeros([nNS,nEW])
	elon = np.zeros([nNS+1,nEW+1])
	elat = np.zeros([nNS+1,nEW+1])		
	for i in np.arange(nNS+1):
		for j in np.arange(nEW+1):
			if i == nNS and j < nEW:
				elon[i,j], elat[i,j] = m(cX[j]-0.5*dX[j], cY[nNS-1]+0.5*dY[nNS-1], inverse=True)
			elif j == nEW and i < nNS:
				elon[i,j], elat[i,j] = m(cX[nEW-1]+0.5*dX[nEW-1], cY[i]-0.5*dY[i], inverse=True)
			elif j == nEW and i == nNS:
				elon[i,j], elat[i,j] = m(cX[nEW-1]+0.5*dX[nEW-1], cY[nNS-1]+0.5*dY[nNS-1], inverse=True)
			else:
				elon[i,j], elat[i,j] = m(cX[j]-0.5*dX[j], cY[i]-0.5*dY[i], inverse=True)

	sizeNS = np.ma.masked_all([nNS,nEW])
	sizeEW = np.ma.masked_all([nNS,nEW])
	for i in np.arange(nNS):
		for j in np.arange(nEW):			
			clon[i,j], clat[i,j] = m(cX[j], cY[i], inverse=True)
			sizeNS[i,j] = dY[i]*1e-3  #km
			sizeEW[i,j] = dX[j]*1e-3  #km
	
	sizeNS[np.where(clon == 1e30)] = np.ma.masked
	sizeNS[np.where(clat == 1e30)] = np.ma.masked
	sizeEW[np.where(clon == 1e30)] = np.ma.masked
	sizeEW[np.where(clat == 1e30)] = np.ma.masked
	
	# Discrete colorbar	
	cmap = cm.jet
	# extract all colors from the .jet map
	cmaplist = [cmap(i) for i in range(cmap.N)]
	# create the new map
	cmap = cmap.from_list('Custom cmap', cmaplist, cmap.N)

	# define the bins and normalize
	bounds = np.arange(4,10,0.5)
	norm = colors.BoundaryNorm(bounds, cmap.N)
	cmap.set_bad(color='w',alpha=0)

	# plot NS pixel size
	m = set_basemap(projection,clon,clat,elon,elat,lon_0=128.2,satellite_height=35785831.0)
	fig = plt.figure()
	im = map_(fig,m,sizeNS,cmap,bounds.min(),bounds.max(),'NS Pixel Size [km]',elon=elon,elat=elat,norm=norm)	

	ax = plt.gca()
	#[left, bottom, width, height]
	cbaxes = fig.add_axes([0.91, 0.147, 0.02, 0.73]) 
	fig.colorbar(im,ax=ax, cax = cbaxes) 
	plt.savefig('gems_NS_pixelsize.png', transparent='true')
	plt.close(fig)

	# plot EW pixel size
	fig = plt.figure()
	bounds = np.arange(9,26)
	norm = colors.BoundaryNorm(bounds, cmap.N)
	im = map_(fig,m,sizeEW,cmap,bounds.min(),bounds.max(),'EW Pixel Size [km]',elon=elon,elat=elat,norm=norm)	

	ax = plt.gca()
	#[left, bottom, width, height]
	cbaxes = fig.add_axes([0.91, 0.147, 0.02, 0.73]) 
	fig.colorbar(im,ax=ax, cax = cbaxes) 
	plt.savefig('gems_EW_pixelsize.png', transparent='true')
	plt.close(fig)  