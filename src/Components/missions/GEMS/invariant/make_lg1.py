#!/usr/bin/env python
"""
Make the lg1 file for GEMS
"""
import numpy as np
from numpy import radians, degrees
from math import radians, cos, sin, asin, sqrt, tan, atan
from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
from matplotlib import cm, colors
from geopy import distance, point

def haversine(lat1, lon1, lat2, lon2):
		"""
		Calculate the great circle distance between two points 
		on the earth (specified in decimal degrees)
		"""
		# convert decimal degrees to radians 
		lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])

		# haversine formula 
		dlon = lon2 - lon1 
		dlat = lat2 - lat1 
		a = sin(dlat/2)**2 + cos(lat1) * cos(lat2) * sin(dlon/2)**2
		c = 2 * asin(sqrt(a)) 
		r = 6371 # Radius of earth in kilometers. Use 3956 for miles
		return c * r


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

		m.drawparallels(np.arange(0.,40.,10.)) # draw parallels
		m.drawmeridians(np.arange(0.,420.,10.)) # draw meridians    

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
	rsphere          = (6378137.00,6356752.3142)
	centerlon = 114.1  #degrees E
	centerlat = 17.04  #degrees N
	fov       = 70  #degree of view East-West
	latmax = 45
	latmin = -5
	lonmin = 75
	lonmax = 145
	nEW   = 700   
	nNS   = 2000
	doVincenty = False

	#-------
	##  End of User Input
	#-------
	m =  Basemap(projection=projection,lon_0=lon_0,resolution=None,
			 rsphere=(6378137.00,6356752.3142),
			 satellite_height = satellite_height)

	Xmin, Ymin = m(lonmin,latmin)
	Xmax, Ymin = m(lonmax,latmin)
	DX = Xmax - Xmin
	dx = DX/nEW

	Xmax, Ymax = m(lonmax,latmax)
	Xmax, Ymin = m(lonmax,latmin)

	DY = Ymax - Ymin
	dy = DY/nNS

	xes = np.linspace(Xmin, Xmax, nEW+1)
	yes = np.linspace(Ymin, Ymax, nNS+1)
	dx  = xes[1] - xes[0]
	dy  = yes[1] - yes[0]
	xes = xes[0:-1] + 0.5*dx
	yes = yes[0:-1] + 0.5*dy

	# Get lat/lon corner points in map projection coordinates
	clon = np.ma.masked_all([nNS,nEW])
	clat = np.ma.masked_all([nNS,nEW])
	elon = np.ma.masked_all([nNS+1,nEW+1])
	elat = np.ma.masked_all([nNS+1,nEW+1])	

	for i in np.arange(nNS):
		for j in np.arange(nEW):			
			#clon[i,j], clat[i,j] = m(Xmin + j*dx + 0.5*dx, Ymin + i*dy + 0.5*dy, inverse=True)
			clon[i,j], clat[i,j] = m(xes[j], yes[i], inverse=True)

	xes = np.linspace(Xmin, Xmax, nEW+1)
	yes = np.linspace(Ymin, Ymax, nNS+1)
	for i in np.arange(nNS+1):  
		for j in np.arange(nEW+1): 
			#elon[i,j], elat[i,j] = m(Xmin + j*dx, Ymin + i*dy, inverse=True)
			elon[i,j], elat[i,j] = m(xes[j], yes[i], inverse=True)

	pixel_top     = np.ma.masked_all([nNS,nEW])
	pixel_bottom  = np.ma.masked_all([nNS,nEW])
	pixel_left    = np.ma.masked_all([nNS,nEW])
	pixel_right   = np.ma.masked_all([nNS,nEW])
	for i in np.arange(nNS):
		for j in np.arange(nEW):
			if doVincenty:
				if not any(s == 1e30 for s in (elat[i+1,j],elon[i+1,j], elat[i+1,j+1],elon[i+1,j+1])):
					pixel_top[i,j]      = distance.distance((elat[i+1,j],elon[i+1,j]), (elat[i+1,j+1],elon[i+1,j+1])).kilometers

				if not any(s == 1e30 for s in (elat[i,j],elon[i,j]    , elat[i,j+1],elon[i,j+1])):
					pixel_bottom[i,j]   = distance.distance((elat[i,j],elon[i,j])    , (elat[i,j+1],elon[i,j+1])).kilometers

				if not any(s == 1e30 for s in (elat[i+1,j],elon[i+1,j], elat[i,j],elon[i,j])):
					pixel_left[i,j]     = distance.distance((elat[i+1,j],elon[i+1,j]), (elat[i,j],elon[i,j])).kilometers

				if not any(s == 1e30 for s in (elat[i,j+1],elon[i,j+1], elat[i+1,j+1],elon[i+1,j+1])):
					pixel_right[i,j]    = distance.distance((elat[i,j+1],elon[i,j+1]), (elat[i+1,j+1],elon[i+1,j+1])).kilometers
			else:
				#haversine(lat1, lon1, lat2, lon2)
				pixel_top[i,j]      = haversine(elat[i+1,j],elon[i+1,j], elat[i+1,j+1],elon[i+1,j+1])
				pixel_bottom[i,j]   = haversine(elat[i,j],elon[i,j]    , elat[i,j+1],elon[i,j+1])
				pixel_left[i,j]     = haversine(elat[i+1,j],elon[i+1,j], elat[i,j],elon[i,j])
				pixel_right[i,j]    = haversine(elat[i,j+1],elon[i,j+1], elat[i+1,j+1],elon[i+1,j+1])


	pixel_top[np.where(clon == 1e30)] = np.ma.masked
	pixel_top[np.where(clat == 1e30)] = np.ma.masked
	pixel_bottom[np.where(clon == 1e30)] = np.ma.masked
	pixel_bottom[np.where(clat == 1e30)] = np.ma.masked	

	pixel_left[np.where(clon == 1e30)] = np.ma.masked
	pixel_left[np.where(clat == 1e30)] = np.ma.masked
	pixel_right[np.where(clon == 1e30)] = np.ma.masked
	pixel_right[np.where(clat == 1e30)] = np.ma.masked

	clon.mask[clon == 1e30] = True
	clat.mask[clat == 1e30] = True
	elon.mask[elat == 1e30] = True
	elat.mask[elat == 1e30] = True

	# Save to netcdf file
	ncOut 								 = Dataset('gems.lg1.invariant.nc4','w',format='NETCDF4_CLASSIC')
	ncOut.institution 		 = 'NASA/Goddard Space Flight Center'
	ncOut.source      		 = 'Global Model and Assimilation Office'
	ncOut.history     		 = 'Created from make_lg1.py domain information from Jhoon Kim'
	ncOut.references  		 = 'none'
	ncOut.comment     		 = "This file contains GEMS geolocation information at native resolution (before coadding)"
	ncOut.contact          = "Patricia Castellanos <patricia.castellanos@nasa.gov>"
	ncOut.conventions 		 = 'CF'
	ncOut.sat_lat     		 = 0.0
	ncOut.sat_lon     		 = lon_0
	ncOut.sat_alt     		 = satellite_height*1e-3
	ncOut.Earth_radius		 = rsphere[0]*1e-3
	ncOut.Earth_flattening = (rsphere[0] - rsphere[1])/rsphere[0]
	ncOut.NCO        			 = '0.2.1'

	ns   = ncOut.createDimension('ns',nNS)
	ew   = ncOut.createDimension('ew',nEW)
	ns_e = ncOut.createDimension('ns_e',nNS+1)
	ew_e = ncOut.createDimension('ew_e',nEW+1)
	pc   = ncOut.createDimension('pix_corner',4)
	time = ncOut.createDimension('time',1)


	# clat
	varobj = ncOut.createVariable('clat','f4',('ns','ew',))
	varobj.long_name = 'pixel center latitude'
	varobj.units     = 'degrees_north'
	varobj.missing_value = 1e15
	varobj[:] = clat

	#clon
	varobj = ncOut.createVariable('clon','f4',('ns','ew',))
	varobj.long_name = 'pixel center longitude'
	varobj.units     = 'degrees_east'
	varobj.missing_value = 1e15
	varobj[:] = clon

	# elat
	varobj = ncOut.createVariable('elat','f4',('ns_e','ew_e',))
	varobj.long_name = 'latitude at pixel edge'
	varobj.units     = 'degrees_north'
	varobj.missing_value = 1e15
	varobj[:] = elat

	# elon
	varobj = ncOut.createVariable('elon','f4',('ns_e','ew_e',))
	varobj.long_name = 'longitude at pixel edge'
	varobj.units     = 'degrees_east'
	varobj.missing_value = 1e15
	varobj[:] = elon

	# GRADS ew
	varobj = ncOut.createVariable('ew','f4',('ew',))
	varobj.long_name = 'pseudo longitude'
	varobj.units     = 'degrees_east'
	#varobj[:] = np.linspace(lonmin,lonmax,nEW)	
	varobj[:] = clon[0.5*nNS,:]

	# GRADS ew
	varobj = ncOut.createVariable('ns','f4',('ns',))
	varobj.long_name = 'pseudo latitude'
	varobj.units     = 'degrees_north'
	#varobj[:] = np.linspace(latmin,latmax,nNS)
	varobj[:] = clat[:,0.5*nEW]


	# pixel size
	varobj = ncOut.createVariable('pix_size','f4',('pix_corner','ns','ew',))
	varobj.long_name = 'pixel size'
	varobj.units     = 'km'
	varobj[0,:,:] = pixel_top
	varobj[1,:,:] = pixel_right
	varobj[2,:,:] = pixel_bottom
	varobj[3,:,:] = pixel_left

	# time
	varobj = ncOut.createVariable('time','f4',('time',))
	varobj.long_name = 'time'
	varobj.units     = 'hours since 2006-06-01 12:00:00'
	varobj.time_increment = int(10000)
	varobj.begin_date     = int(20060601)
	varobj.begin_time     = int(120000)

	ncOut.close()

	# Discrete colorbar	
	cmap = cm.jet
	# extract all colors from the .jet map
	cmaplist = [cmap(i) for i in range(cmap.N)]
	# create the new map
	cmap = cmap.from_list('Custom cmap', cmaplist, cmap.N)

	# --
	# EW
	# --
	# define the bins and normalize
	bounds = np.arange(7,10.5,0.5)
	norm = colors.BoundaryNorm(bounds, cmap.N)
	cmap.set_bad(color='w',alpha=0)

	m = set_basemap(projection,clon,clat,elon,elat,lon_0=128.2,satellite_height=35785831.0)
	fig = plt.figure()
	fig.set_size_inches(18.5, 10.5)
	im = map_(fig,m,pixel_bottom,cmap,bounds.min(),bounds.max(),'EW Pixel Size [km] nX = {}'.format(nEW),elon=elon,elat=elat,norm=norm)
	
	ax = plt.gca()
	#[left, bottom, width, height]
	cbaxes = fig.add_axes([0.91, 0.147, 0.02, 0.73]) 
	fig.colorbar(im,ax=ax, cax = cbaxes) 
	plt.savefig('gems_EW_pixelsize_nX={}.png'.format(nEW), transparent='true')
	plt.close(fig)  

	# --
	# NS
	# --
	# define the bins and normalize
	bounds = np.arange(2,7.5,0.5)
	norm = colors.BoundaryNorm(bounds, cmap.N)
	cmap.set_bad(color='w',alpha=0)
	fig = plt.figure()
	fig.set_size_inches(18.5, 10.5)
	im = map_(fig,m,pixel_right,cmap,bounds.min(),bounds.max(),'NS Pixel Size [km] nY = {}'.format(nNS),elon=elon,elat=elat,norm=norm)
	
	ax = plt.gca()
	#[left, bottom, width, height]
	cbaxes = fig.add_axes([0.91, 0.147, 0.02, 0.73]) 
	fig.colorbar(im,ax=ax, cax = cbaxes) 
	plt.savefig('gems_NS_pixelsize_nY={}.png'.format(nNS), transparent='true')
	plt.close(fig)  	

	# # --
	# # Area
	# # --
	# # define the bins and normalize
	# bounds = np.arange(10,90,10)
	# norm = colors.BoundaryNorm(bounds, cmap.N)
	# cmap.set_bad(color='w',alpha=0)
	# fig = plt.figure()
	# fig.set_size_inches(18.5, 10.5)
	# im = map_(fig,m,pixel_left*pixel_top,cmap,bounds.min(),bounds.max(),'Pixel Area [km^2]',elon=elon,elat=elat,norm=norm)
	
	# ax = plt.gca()
	# #[left, bottom, width, height]
	# cbaxes = fig.add_axes([0.91, 0.147, 0.02, 0.73]) 
	# fig.colorbar(im,ax=ax, cax = cbaxes) 
	# plt.savefig('gems_pixelarea.png', transparent='true')
	# plt.close(fig)  		