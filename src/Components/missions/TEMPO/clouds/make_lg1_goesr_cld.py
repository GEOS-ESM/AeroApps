#!/usr/bin/env python3
"""
Make the lg1 file for GOES-R type resolution overlaid on TEMPO with GEOS-R parking spot
"""
import numpy as np
from numpy import radians, degrees
from math import radians, cos, sin, asin, sqrt, tan, atan
from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
from matplotlib import cm, colors
from geopy import distance, point
from dateutil.parser import parse         as isoparser
from datetime        import timedelta
def haversine(lat1, lon1, lat2, lon2):
        """
        Calculate the great circle distance between two points 
        on the earth (specified in decimal degrees)
        """
        # convert decimal degrees to radians 
        lon1, lat1, lon2, lat2 = list(map(radians, [lon1, lat1, lon2, lat2]))

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

        meridian  = np.arange(-125.,60.,5.)
        parallels = np.arange(-10.,80.,5.)
        m.drawparallels(parallels) # draw parallels
        m.drawmeridians(meridian) # draw meridians    

        meridian  = np.delete(meridian,len(meridian)*0.5)
        for i in np.arange(len(meridian)):
            plt.annotate(np.str(meridian[i]),xy=m(meridian[i],20),xycoords='data')
        for i in np.arange(len(parallels)):
            plt.annotate(np.str(parallels[i]),xy=m(110,parallels[i]),xycoords='data')
  

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

        print()
        print(name)
        print('Native     Bounding box: ', BBOX) 
        print('Recentered Bounding box: ', Bbox) 
        print('Normalized Bounding box: ', bbox) 

        return (X,Y,BBOX,Bbox)

#---
def scanTimes(nNS,nEW,tBeg):
    """
    Calculate TEMPO scantimes. 
    """
    tScan = np.array([tBeg + i * timedelta(seconds=2.85) for i in np.arange(nEW) ])
    tScan = tScan[-1::-1] # flip times
    return np.tile(tScan,(nNS,1))       

#---
def writeNC(elon,elat,clon,clat,pixel_top,pixel_bottom,pixel_right,pixel_left,DT,fname):
    # Save to netcdf file
    ncOut                                = Dataset(fname,'w',format='NETCDF4_CLASSIC')
    ncOut.institution        = 'NASA/Goddard Space Flight Center'
    ncOut.source             = 'Global Model and Assimilation Office'
    ncOut.history            = 'Created from make_lg1_tempo_cld.py'
    ncOut.references         = 'none'
    ncOut.comment            = "This file contains TEMPO geolocation information for cloud modeling"
    ncOut.contact          = "Patricia Castellanos <patricia.castellanos@nasa.gov>"
    ncOut.conventions        = 'CF'
    ncOut.sat_lat            = 0.0
    ncOut.sat_lon            = lon_0
    ncOut.sat_alt            = satellite_height*1e-3
    ncOut.Earth_radius       = rsphere[0]*1e-3
    ncOut.Earth_flattening = (rsphere[0] - rsphere[1])/rsphere[0]
    ncOut.NCO                    = '0.2.1'


    nNS,nEW = clon.shape
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
    varobj.missing_value = 1e15
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

    # scanTime
    varobj = ncOut.createVariable('scanTime','f4',('ew',))
    varobj.units = "seconds since 0001-01-01 00:00:00.0"
    varobj[:] = DT

    ncOut.close()





#---
if __name__ == '__main__':

    outFile    = 'goes-r.lg1.cld.invariant.'
    inFile     = 'tempo.lg1.cld.invariant.'
    layout     = '41'
    projection = 'geos'
    lon_0      = -75
    lat_0      = 0.0


    satellite_height = 35785831.0
    rsphere          = (6378137.00,6356752.3142)

    #-------
    ##  End of User Input
    #-------
    m =  Basemap(projection=projection,lon_0=lon_0,resolution=None,
             rsphere=(6378137.00,6356752.3142),
             satellite_height = satellite_height)

    # Read in TEMPO
    ntiles = int(layout[0])*int(layout[1])
    for tile in range(ntiles):
        print('Reading file',inFile + layout + str(tile) + '.nc4')
        ncTempo = Dataset(inFile + layout + str(tile) + '.nc4')
        clon = ncTempo.variables['clon'][:]
        clat = ncTempo.variables['clat'][:]
        pixel_top    = np.squeeze(ncTempo.variables['pix_size'][0,:,:])
        pixel_right  = np.squeeze(ncTempo.variables['pix_size'][1,:,:])
        pixel_bottom = np.squeeze(ncTempo.variables['pix_size'][2,:,:])
        pixel_left   = np.squeeze(ncTempo.variables['pix_size'][3,:,:])

        if not hasattr(clon,'mask'):
            print('tile',tile)
            clon = np.ma.array(clon,mask=np.zeros(clon.shape).astype(bool))

        if not hasattr(clat,'mask'):
            clat = np.ma.array(clat,mask=np.zeros(clat.shape).astype(bool))

        if not hasattr(pixel_top,'mask'):
            pixel_top = np.ma.array(pixel_top,mask=np.zeros(pixel_top.shape).astype(bool))

        if not hasattr(pixel_right,'mask'):
            pixel_right = np.ma.array(pixel_right,mask=np.zeros(pixel_right.shape).astype(bool))

        if not hasattr(pixel_bottom,'mask'):
            pixel_bottom = np.ma.array(pixel_bottom,mask=np.zeros(pixel_bottom.shape).astype(bool))

        if not hasattr(pixel_left,'mask'):
            pixel_left = np.ma.array(pixel_left,mask=np.zeros(pixel_left.shape).astype(bool))

        if tile == 0:
            PP_bottom = np.ma.array(pixel_bottom,mask=np.zeros(pixel_bottom.shape).astype(bool))
        else:
            PP_bottom = np.ma.append(PP_bottom,pixel_bottom,axis=1)


        X, Y = m(clon,clat)
        I = (X == 1e30) | (Y == 1e30)
        X.mask[I] = True
        Y.mask[I] = True
        clon.mask = X.mask | Y.mask
        clat.mask = X.mask | Y.mask
        pixel_top.mask    = X.mask | Y.mask
        pixel_bottom.mask = X.mask | Y.mask
        pixel_left.mask   = X.mask | Y.mask
        pixel_right.mask  = X.mask | Y.mask

        elon = ncTempo.variables['elon'][:]
        elat = ncTempo.variables['elat'][:]

        if not hasattr(elon,'mask'):
            elon = np.ma.array(elon,mask=np.zeros(elon.shape).astype(bool))

        if not hasattr(elat,'mask'):
            elat = np.ma.array(elat,mask=np.zeros(elat.shape).astype(bool))

        if tile == 0:
            EElon     = np.ma.array(elon,mask=np.zeros(elon.shape).astype(bool))
        else:
            EElon     = np.ma.append(EElon[:,0:-1],elon,axis=1)
    
        X, Y = m(elon,elat)
        I = (X == 1e30) | (Y == 1e30)
        X.mask[I] = True
        Y.mask[I] = True        
        elon.mask = X.mask | Y.mask
        elat.mask = X.mask | Y.mask

        mask_1 = elon.mask[0:-1,0:-1]
        mask_2 = elon.mask[0:-1,1:]
        mask_3 = elon.mask[1:,0:-1]
        mask_4 = elon.mask[1:,1:]

        clon.mask         = mask_1 | mask_2 | mask_3 | mask_4
        clat.mask         = mask_1 | mask_2 | mask_3 | mask_4
        pixel_top.mask    = mask_1 | mask_2 | mask_3 | mask_4
        pixel_bottom.mask = mask_1 | mask_2 | mask_3 | mask_4
        pixel_left.mask   = mask_1 | mask_2 | mask_3 | mask_4
        pixel_right.mask  = mask_1 | mask_2 | mask_3 | mask_4



        DT = ncTempo.variables['scanTime'][:]
        print('writing file',outFile + layout + str(tile) + '.nc4')
        writeNC (elon, elat, clon, clat, pixel_top, pixel_bottom, pixel_right, pixel_left, DT,
                  fname = outFile + layout + str(tile) + '.nc4')

        if tile == 0:
            P_bottom = pixel_bottom
            P_right  = pixel_right
            Elon     = elon
            Elat     = elat
        else:
            P_bottom = np.ma.append(P_bottom,pixel_bottom,axis=1)
            P_right  = np.ma.append(P_right,pixel_right,axis=1)
            Elon     = np.ma.append(Elon[:,0:-1],elon,axis=1)
            Elat     = np.ma.append(Elat[:,0:-1],elat,axis=1)

        ncTempo.close()


    # ---
    # Make some plots
    # ---
    pixel_bottom = P_bottom
    pixel_right  = P_right
    elon         = Elon
    elat         = Elat
    lon_0        = -100
    nNS = pixel_bottom.shape[0]
    nEW = pixel_bottom.shape[1]

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
    bounds = np.arange(0,10,1)
    norm = colors.BoundaryNorm(bounds, cmap.N)
    cmap.set_bad(color='w',alpha=0)

    m = set_basemap(projection,clon,clat,elon,elat,lon_0=lon_0,satellite_height=35785831.0)
    fig = plt.figure()
    fig.set_size_inches(18.5, 10.5)
    im = map_(fig,m,pixel_bottom,cmap,bounds.min(),bounds.max(),'EW Pixel Size [km] nX = {}'.format(nEW),elon=elon,elat=elat,norm=norm)
    
    ax = plt.gca()
    #[left, bottom, width, height]
    cbaxes = fig.add_axes([0.91, 0.147, 0.02, 0.73]) 
    fig.colorbar(im,ax=ax, cax = cbaxes) 
    plt.savefig('goes-r_EW_pixelsize_nX={}.png'.format(nEW), transparent='true')
    plt.close(fig)  

    # --
    # NS
    # --
    # define the bins and normalize
    bounds = np.arange(0,5.5,0.5)
    norm = colors.BoundaryNorm(bounds, cmap.N)
    cmap.set_bad(color='w',alpha=0)
    fig = plt.figure()
    fig.set_size_inches(18.5, 10.5)
    im = map_(fig,m,pixel_right,cmap,bounds.min(),bounds.max(),'NS Pixel Size [km] nY = {}'.format(nNS),elon=elon,elat=elat,norm=norm)
    
    ax = plt.gca()
    #[left, bottom, width, height]
    cbaxes = fig.add_axes([0.91, 0.147, 0.02, 0.73]) 
    fig.colorbar(im,ax=ax, cax = cbaxes) 
    plt.savefig('goes-r_NS_pixelsize_nY={}.png'.format(nNS), transparent='true')
    plt.close(fig)      
