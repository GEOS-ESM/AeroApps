#!/usr/bin/env python3

"""
    Utility to regrid and gap fill MCD43C1 Collections onto a geostationary domain.

    P. Castellanos, Nov 2015
"""

import os
from optparse        import OptionParser
from datetime        import datetime, timedelta
from dateutil.parser import parse         as isoparser
import numpy as np
import mcd43c
from netCDF4 import Dataset
from glob     import glob
from collections import defaultdict
import calendar

MISSING = -99999
MODIS_MISSING = 32.767
Kernels = 3
kernel_names = ['KISO','KVOL','KGEO']
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

outbands = ('b1_645',
            'b2_856',
            'b3_465',
            'b4_553',
            'b5_1241',
            'b6_1629',
            'b7_2114')

bandname = {}
for i, b in enumerate(outbands):
  bandname[b] = 'Band'+str(i+1)

#---
def getCoords(geoFile,verbose):
    """
    Get coordinates from geo-location file.
    Assumes the file has albed been tighten in the E-W domain
    """
    if verbose:
       print(" <> Getting GEO coordinates from ", geoFile)
    nc = Dataset(geoFile)
    lon = nc.variables['clon'][:,:]
    lat = nc.variables['clat'][:,:]
    missing = nc.variables['clon'].missing_value
    return (nc,lon,lat,missing)

#----
def _copyVar(ncIn,ncOut,name,dtype='f4',layout=None,zlib=False):
    """
    Create variable *name* in output file and copy its
    content over,
    """
    x = ncIn.variables[name]
    print(x.dimensions)
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
        print('_copyVar',nX,nY,i,i_,j_)


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
        raise ValueError("invalid rank of <%s>: %d"%(name,rank))


#----
def writenc(mcdData,ncGeo,clon,clat,options,xchunk=150, ychunk=200):
    """
    Write a NetCDF file with sampled MCD43C1 kernel weights on geostationary grid
    described by (clon,clat).
    """

    # Gridded Dimensions
    # ------------------
    nNS, nEW = clon.shape

    # Open NC file
    # ------------
    nc = Dataset(options.outFile,'w',format=options.format)

    # Set global attributes
    # ---------------------
    nc.title = options.title
    nc.institution = 'NASA/Goddard Space Flight Center'
    nc.source = 'Global Model and Assimilation Office'
    nc.history = 'Created from MCD43C1 v005 collections by mcd43c_sampler.py'
    nc.references = 'n/a'
    nc.comment = 'This file contains BRDF Kernels weights for the RTLS model for 8 MODIS bands sampled on a geostationary grid'
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
    k = nc.createDimension('Kernels',Kernels)

    # Add pseudo dimensions for GrADS compatibility
    # -------------------------------------------
    _copyVar(ncGeo,nc,'ew',dtype='f4',zlib=False)
    _copyVar(ncGeo,nc,'ns',dtype='f4',zlib=False)

    # Save lon/lat if so desired
    # --------------------------
    if options.coords:
        _copyVar(ncGeo,nc,'clon',dtype='f4',zlib=False)
        _copyVar(ncGeo,nc,'clat',dtype='f4',zlib=False)

    # Loop over Bands writing each dataset
    #---------------------------------------
    dim = ('Kernels','ns','ew')
    chunks = (1,ychunk, xchunk)
    for b in outbands:
      this = nc.createVariable(bandname[b],'f4',dim,
                                     zlib=options.zlib,
                                     chunksizes=chunks)  

      this.long_name = bandname[b] + ' BRDF Kernel weight: isotropic, volumetric, geometric'
      this.missing_value = -99999
      this.unit = 'none'  

      data  = np.ma.masked_all([Kernels,nNS,nEW])
      temp = np.ma.masked_all([nNS,nEW])    
      for i,k in enumerate(kernel_names):      
        temp[~clon.mask] = getattr(mcdData,k + '_' + b)    
        data[i,:,:] = temp

      this[:] = data


    nc.close()  


#----
def season_dic(years):

  if type(years) is int:
    years = [years]

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

  return dates

#----
def gap_fill(mcdData, clon, clat, FRLAND, season, ncClim):
  """
    Gap fill BRDF kernel data with seasonal average. 
  """
  nNS, nEW = clon.shape 
  for b in bands:
    bdata = getattr(mcdData,b)
    sdata = ncClim.variables[season + '_' + b][:]

    data  = np.ma.masked_all([nNS,nEW])
    data[~clon.mask] = bdata    
    data[data > 32]    = sdata[data > 32]
    data[FRLAND< 0.99] = np.ma.masked

    setattr(mcdData,b,data[~clon.mask])



#------------------------------------ M A I N ------------------------------------

if __name__ == "__main__":

  inst     = 'SENTINEL-4'
  title    = 'MCD43C1 SENTINEL-4 Sampler'
  format   = 'NETCDF4_CLASSIC'
  outFile  = 'sentinel-4.mcd43c_sampler.nc4'
  landFile = 'sentinel-4-g5nr.lb2.asm_Nx.nc4'
  climFile = 'sentinel-4.mcd43c1_climatology.nc4'

  # Geolocation default
  # -------------------
  calculon = '/nobackup/'+inst+'/LevelG/invariant/'
  nccs = '/discover/nobackup/projects/gmao/osse2/pub/c1440_NR/OBS/'+inst+'/DATA/LevelG/invariant/'
  if os.path.exists(calculon): 
      geoFile      = calculon + inst.lower()+'.lg1.invariant.nc4'
      datadir      = '/nobackup/3/pcastell/MODIS/MCD43C1/'
      landFile     = '/nobackup/'+inst+'/LevelB/invariant/'+inst.lower()+'-g5nr.lb2.asm_Nx.nc4'
  elif os.path.exists(nccs): 
      geoFile      = nccs  +inst.lower()+'.lg1.invariant.nc4'
  else:
      geoFile      = inst.lower()+'.lg1.invariant.nc4'
      datadir      = '/nobackup/3/pcastell/MODIS/MCD43C1/'


#   Parse command line options
#   --------------------------
  parser = OptionParser(usage="Usage: %prog [OPTIONS] t1 [t2 t3...] ",
                        version='1.0.0' )

  parser.add_option("-I", "--isotime",
                    action="store_true", dest="iso",
                    help="Input is isotime. Default is YYYYJJJ.")  

  parser.add_option("-C", "--coords",
                    action="store_true", dest="coords", default=True,
                    help="Save 2D lon/lat coordinates on output file.")    

  parser.add_option("-g", "--geolocation", dest="geoFile", default=geoFile,
            help='Level G1 geo-location file (default=%s)'%geoFile)

  parser.add_option("-o", "--output", dest="outFile", default=outFile,
            help="Output NetCDF file (default=%s)"\
                        %outFile )

  parser.add_option("-l", "--landfile", dest="landFile", default=landFile,
            help="Land Fraction NetCDF file (default=%s)"\
                        %landFile )    

  parser.add_option("-t", "--title", dest="title", default=title,
            help="Output file title, typically the collection name (default=%s)"\
                        %title )

  parser.add_option("-v", "--verbose",
                    action="store_true", dest="verbose",
                    help="Verbose mode.")

  parser.add_option("-i", "--instrument",
                dest="inst", default=inst, help="Instrument name (default=%s)"\
                        %inst)                                            

  parser.add_option("-f", "--format", dest="format", default=format,
            help="Output file format: one of NETCDF4, NETCDF4_CLASSIC, NETCDF3_CLASSIC or NETCDF3_64BIT (default=%s)"%format )

  parser.add_option("-N", "--nozip",
                    action="store_true", dest="nozip",
                    help="Do not shave/compress output files.")


  (options, args) = parser.parse_args()

  t1 = []
  if options.iso:
      for t in args:
          t1.append(isoparser(t))
  else:
      for t in args:
          t1.append(datetime.strptime(t, '%Y%j'))

  # Create consistent file name extension
  # -------------------------------------
  name, ext = os.path.splitext(options.outFile)
  if 'NETCDF4' in options.format:
      options.outFile = name + '.nc4'
  elif 'NETCDF3' in options.format:
      options.outFile = name + '.nc'
  else:
      raise ValueError('invalid extension <%s>'%ext)
  options.zlib = not options.nozip

  # Create (x,y,t) coordinates
  # --------------------------
  ncGeo, clon, clat, MISSING = getCoords(options.geoFile,options.verbose)       

  # Get Land Fraction
  #---------------------
  ncLand = Dataset(options.landFile)
  FRLAND = np.squeeze(ncLand.variables['FRLAND'][:])
  ncLand.close()  

  ncClim = Dataset(climFile)     
 
  for t in t1:
    dates = season_dic(t.year)
    path = glob(datadir + '*' + t.strftime('%Y%j') + '*')
    for fn in path: assert os.path.exists(fn), fn + ' DOES NOT EXIST' 
    mcdData = mcd43c.McD43C(path,clon[~clon.mask],clat[~clat.mask],Verb=0)

    # Gap Fill Data
    #----------------
    for s in dates:
      if t in dates[s]:
        season = s

    gap_fill(mcdData,clon,clat,FRLAND,season,ncClim)

    writenc(mcdData,ncGeo,clon,clat,options)

  #ncClim.close()
  #ncGeo.close()

"""
elon = ncGeo.variables[u'elon'][:,:]
elat = ncGeo.variables[u'elat'][:,:]
projection = 'geos'
m = set_basemap(projection,clon,clat,elon,elat)
a = mcdData.KISO_b1_645
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
