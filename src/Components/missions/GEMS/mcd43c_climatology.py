#!/usr/bin/env python

"""
    Utility to write a global MCD43C1 seasonal climatology.
    Will be used to gapfill 16 day MCD43C1 files

    P. Castellanos, Nov 2015
"""

import os
from optparse        import OptionParser
from datetime        import datetime, timedelta
import calendar
from dateutil.parser import parse         as isoparser
import numpy as np
import mcd43c
from netCDF4 import Dataset
from glob     import glob
from collections import defaultdict

MISSING = 32.767
Kernels = 3
bands = ('KISO_b1_645',
         'KISO_b2_856',
         'KISO_b3_465',
         'KISO_b4_553',
         'KISO_b5_1241',
         'KISO_b6_1629',
         'KISO_b7_2114',
         'KVOL_b1_645',
         'KVOL_b2_856',
         'KVOL_b3_465',
         'KVOL_b4_553',
         'KVOL_b5_1241',
         'KVOL_b6_1629',
         'KVOL_b7_2114',
         'KGEO_b1_645',
         'KGEO_b2_856',
         'KGEO_b3_465',
         'KGEO_b4_553',
         'KGEO_b5_1241',
         'KGEO_b6_1629',
         'KGEO_b7_2114')
seasons = 'DJF', 'MAM', 'JJA', 'SON'
#----
def _copyVar(ncIn,ncOut,name,dtype='f4',layout=None,zlib=False):
    """
    Create variable *name* in output file and copy its
    content over,
    """
    x = ncIn.variables[name]
    print x.dimensions
    y = ncOut.createVariable(name,dtype,x.dimensions,zlib=zlib)
    y.long_name = x.long_name
    y.units = x.units 
    try:
        y.missing_value = x.missing_value
    except:
        pass
    rank = len(x.shape)

    if layout is not None:
        nX = int(layout[0])
        nY = int(layout[1])
        i  = int(layout[2])

        i_ = i%nX
        j_ = int(i/nX) 
        print '_copyVar',nX,nY,i,i_,j_


    if rank == 1:
        if layout is None:
            y[:] = x[:]
        else:
            nEW, = x.shape
            Xstart = i_*nEW/nX
            Xend   = Xstart + nEW/nX
            y[:] = x[Xstart:Xend]
    elif rank == 2:
        if layout is None:
            y[:,:] = x[:,:]
        else:
            nNS, nEW = x.shape 
            Xstart = i_*nEW/nX
            Xend   = Xstart + nEW/nX
            Ystart = j_*nNS/nY
            Yend   = Ystart + nNS/nY

            y[:,:] = x[Ystart:Yend,Xstart:Xend]
    elif rank == 3:
        if layout is None:
            y[:,:,:] = x[:,:,:]
        else:
            nNS, nEW = x.shape 
            Xstart = i_*nEW/nX
            Xend   = Xstart + nEW/nX
            Ystart = j_*nNS/nY
            Yend   = Ystart + nNS/nY

            y[:,:,:] = x[:,Ystart:Yend,Xstart:Xend]
    else:
        raise ValueError, "invalid rank of <%s>: %d"%(name,rank)
#----
def writenc(Avgs, Cnts, ncGeo, clon, clat, outFile):
    """
    Write a NetCDF file with MCD43C1 seasonal climatology
    [DJF, MAM, JJA, SON]
    """

    # Gridded Dimensions
    # ------------------
    nNS, nEW = clon.shape

    # Open NC file
    # ------------
    nc = Dataset(outFile,'w',format='NETCDF4_CLASSIC')

    # Set global attributes
    # ---------------------
    nc.title = 'MCD43C1 Seasonal Climatology'
    nc.institution = 'NASA/Goddard Space Flight Center'
    nc.source = 'Global Model and Assimilation Office'
    nc.history = 'Created from MCD43C1 v005 collections by mcd43c_climatology.py'
    nc.references = 'n/a'
    nc.comment = 'This file contains seasonal average BRDF Kernels weights for the RTLS model for 8 MODIS bands sampled on a geostationary grid'
    nc.contact = 'Patricia Castellanos <patricia.castellanos@nasa.gov>'
    nc.Conventions = 'CF'    
    nc.BAND1 = "620-670nm"
    nc.BAND2 = "841-875nm"
    nc.BAND3 = "459-479nm"
    nc.BAND4 = "545-565nm"
    nc.BAND5 = "1230-1250nm"
    nc.BAND6 = "1628-1652nm"
    nc.BAND7 = "2105-2155nm"

    # Create dimensions
    # -----------------
    x = nc.createDimension('ew',nEW)
    y = nc.createDimension('ns',nNS)

    # Add pseudo dimensions for GrADS compatibility
    # -------------------------------------------
    _copyVar(ncGeo,nc,u'ew',dtype='f4',zlib=False)
    _copyVar(ncGeo,nc,u'ns',dtype='f4',zlib=False)

    # Save lon/lat 
    # --------------------------
    _copyVar(ncGeo,nc,u'clon',dtype='f4',zlib=False)
    _copyVar(ncGeo,nc,u'clat',dtype='f4',zlib=False)

    # Loop over Bands writing each dataset
    #---------------------------------------
    dim = ('ns','ew',)
    for s in seasons:
      for b in bands:
        this = nc.createVariable(s + '_' + b,'f4',dim)  

        this.long_name = s + ' ' + b + ' Kernel weight'
        this.missing_value = -99999
        this.unit = 'none'  

        data  = np.ma.masked_all([nNS,nEW])
        ava = np.ma.array(getattr(Avgs[s],b))
        cnt = getattr(Cnts[s],b)
        ava[cnt == 0] = np.ma.masked
        data[~clon.mask] = ava

        this[:] = data


    nc.close()  

# -----
def parse_date(filelist):

  filedates = np.empty(0)
  for fn in filelist:
    fn = fn.split('/')[-1]
    fn = fn.split('A')[1]
    fn = fn.split('.')[0]
    filedates = np.append(filedates,datetime.strptime(fn, '%Y%j'))

  return filedates

# -----
def getCoords(geoFile,verbose):
    """
    Get coordinates from geo-location file.
    Assumes the file has albed been tighten in the E-W domain
    """
    if verbose:
       print " <> Getting GEO coordinates from ", geoFile
    nc = Dataset(geoFile)
    lon = nc.variables[u'clon'][:,:]
    lat = nc.variables[u'clat'][:,:]
    missing = nc.variables[u'clon'].missing_value
    return (nc,lon,lat,missing)
#------------------------------------ M A I N ------------------------------------

if __name__ == "__main__":

    format = 'NETCDF4_CLASSIC'
    outFile = 'gems.mcd43c1_climatology.nc4'
    geoFile = 'gems.lg1.invariant.nc4'

    path    = '/nobackup/3/pcastell/MODIS/MCD43C1/'
    years   = 2006,2007


    # Create Dictionaries of dates for each season
    #-----------------------------------------------
    months = {'DJF': (12,1,2), 'MAM': (3,4,5), 'JJA': (6,7,8), 'SON': (9,10,11)}
    dates   = {}
    dates   = defaultdict(lambda: np.empty(0), dates)
    for s in months:
      for y in years:
        for m in months[s]:
          num_days = calendar.monthrange(y, m)[1]
          newdates  = [datetime(y, m, day) for day in range(1, num_days+1)]
          dates[s] = np.append(dates[s],newdates)

    # Create a list of the files and their dates
    #-------------------------------------------
    filelist = glob(path + '*.hdf')
    filedates = parse_date(filelist)

    # Open each file
    # Grid to GEMS and
    # add to appropriate season
    #-------------------------------
    ncGrid, clon, clat, GEOMISSING = getCoords(geoFile,False)

    Avgs   = defaultdict(lambda: mcd43c.McD43C(filelist[0],clon[~clon.mask],clat[~clat.mask],Verb=0))
    Cnts   = defaultdict(lambda: mcd43c.McD43C(filelist[0],clon[~clon.mask],clat[~clat.mask],Verb=0))
    for s in seasons:
      for b in bands:
        setattr(Avgs[s],b,np.zeros(len(clon[~clon.mask])))
        setattr(Cnts[s],b,np.zeros(len(clon[~clon.mask])))

    for index, fn in enumerate(filelist):
      mcdData = mcd43c.McD43C(fn,clon[~clon.mask],clat[~clat.mask],Verb=0)
      for s in seasons:
        if filedates[index] in dates[s]:
          for b in bands:
            data = getattr(mcdData,b)
            Avgs[s].__dict__[b][data < 32]  += data[data < 32]
            Cnts[s].__dict__[b][data < 32]  += 1 

    # Get average; divide by counts
    # ------------------------------
    for s in seasons:
      for b in bands:
        cnt = Cnts[s].__dict__[b]
        Avgs[s].__dict__[b][cnt >0] = Avgs[s].__dict__[b][cnt >0] / cnt[cnt >0]

    writenc(Avgs,Cnts,ncGrid,clon,clat,outFile)

"""
elon = ncGeo.variables[u'elon'][:,:]
elat = ncGeo.variables[u'elat'][:,:]
projection = 'geos'
m = set_basemap(projection,clon,clat,elon,elat)
a = mcdData.BRDF_Albedo_Parameter1_Band1
nNS, nEW = clon.shape
data = np.ma.masked_all([nNS,nEW])
data[~clon.mask] = a
data[data>32] = np.ma.masked
norm = None
minval = 0
maxval = 1
cmap = cm.jet
cmap.set_bad(color='w',alpha=0)
fig = plt.figure()
fig.set_size_inches(18.5, 10.5)
fig.set_dpi(100)
im = m.imshow(data,cmap=cmap,vmin=minval,vmax=maxval,norm=norm)
m.drawcoastlines(color='purple')
m.drawcountries(color='purple')
m.drawstates(color='purple')
plt.savefig('test_mcd43.png')

from matplotlib import cm, colors
cmap = cm.jet
cmap.set_bad(color='w',alpha=0)
lon = np.arange(-180,180,1)
lat = np.arange(-90,90,1)
lon,lat = np.meshgrid(lon,lat)
mcdData = mcd43c.McD43C(path,lon.flatten(),lat.flatten(),Verb=0)
a = mcdData.BRDF_Albedo_Parameter1_Band1
a = np.ma.array(a)
a[a>32] = np.ma.masked
import matplotlib.pyplot as plt
plt.contourf(a.reshape(lon.shape))

from mpl_toolkits.basemap import Basemap
m = Basemap()
m.contourf(lon,lat,a.reshape(lon.shape),latlon=True)
m.drawcoastlines(color='purple')
"""
