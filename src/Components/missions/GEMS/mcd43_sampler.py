#!/usr/bin/env python

"""
    Utility to regrid MCD43B1 Collections onto GEMS domain.

    P. Castellanos, Nov 2015
"""

import os
from optparse        import OptionParser
from datetime        import datetime, timedelta
from dateutil.parser import parse         as isoparser
import numpy as np
from pyobs import mcd43
from netCDF4 import Dataset
from scipy.interpolate import interp2d

MISSING = -99999
MODIS_MISSING = 32.767
Kernels = 3
bands = {'BRDF_b1_645'  : 'Band1',
          'BRDF_b2_856' : 'Band2',
          'BRDF_b3_465' : 'Band3',
          'BRDF_b4_553' : 'Band4',
          'BRDF_b5_1241': 'Band5',
          'BRDF_b6_1629': 'Band6',
          'BRDF_b7_2114': 'Band7'}

#---
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
def writenc(mcdData,ncGeo,clon,clat,options,xchunk=150, ychunk=200):
    """
    Write a NetCDF file with sampled MCD43B1 kernel weights on geostationary grid
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
    nc.history = 'Created from MCD43B1 v005 collections by mcd43_sampler.py'
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
    _copyVar(ncGeo,nc,u'ew',dtype='f4',zlib=False)
    _copyVar(ncGeo,nc,u'ns',dtype='f4',zlib=False)

    # Save lon/lat if so desired
    # --------------------------
    if options.coords:
        _copyVar(ncGeo,nc,u'clon',dtype='f4',zlib=False)
        _copyVar(ncGeo,nc,u'clat',dtype='f4',zlib=False)

    # Loop over Bands writing each dataset
    #---------------------------------------
    dim = ('Kernels','ns','ew')
    chunks = (1,ychunk, xchunk)
    for b in bands:
      this = nc.createVariable(bands[b],'f4',dim,
                                     zlib=options.zlib,
                                     chunksizes=chunks)  

      this.long_name = bands[b] + ' BRDF Kernel weight: isotropic, volumetric, geometric'
      this.missing_value = -99999
      this.unit = 'none'  

      data  = np.ma.masked_all([Kernels,nNS,nEW])
      temp = np.ma.masked_all([nNS,nEW])    
      bdata = getattr(mcdData,b)
      for k in np.arange(Kernels):            
        temp[~clon.mask] = bdata[:,k]
        data[k,:,:] = temp

      this[:] = data


    nc.close()  

#----
def gap_fill(mcdData, clon, clat, options):
  """
    Gap fill BRDF kernel data with linear interpolation. 
  """

  # Get Land Fraction
  #---------------------
  ncLand = Dataset(options.landFile)
  FRLAND = np.squeeze(ncLand.variables[u'FRLAND'][:])
  ncLand.close()


  for b in bands:
    bdata = getattr(mcdData,b)
    kdata = np.ma.masked_all([nNS,nEW])
    for k in np.arange(Kernels):
      kdata[~clon.mask] = bdata[:,k]
      masktrain  = np.where((FRLAND >= 0.99) & (~clon.mask) & (abs(kdata - MISSING) > 0.01) & (abs(kdata - MODIS_MISSING) > 0.01))
      maskinterp = np.where((FRLAND >= 0.99) & (~clon.mask) & (abs(kdata - MISSING) <= 0.01) & (abs(kdata - MODIS_MISSING) <= 0.01))
      f     = interp2d(clon[masktrain], clat[masktrain], kdata[masktrain])

      kdata[maskinterp] = f(clon[maskinterp],clat[maskinterp])
      bdata[:,k]        = kdata[~clon.mask]

    setattr(mcdData,b,bdata)



#------------------------------------ M A I N ------------------------------------

if __name__ == "__main__":

    inst   = 'GEMS'
    title  = 'MCD43B1 GEMS Sampler'
    format = 'NETCDF4_CLASSIC'
    rcFile  = 'geo_sampler_omi_ler.rc'
    outFile = 'gems_mcd43_sampler.nc4'
    landFile = 'gems-g5nr.lb2.asm_Nx.nc4'
    layout       = None   

    # Geolocation default
    # -------------------
    calculon = '/nobackup/'+inst+'/LevelG/invariant/'
    nccs = '/discover/nobackup/projects/gmao/osse2/pub/c1440_NR/OBS/'+inst+'/DATA/LevelG/invariant/'
    if os.path.exists(calculon): 
        geoFile      = calculon + inst.lower()+'.lg1.invariant.nc4'
        scanTimeFile = calculon + inst.lower()+'.lg1.scantime.nc4'
        datadir      = '/nobackup/3/pcastell/MODIS/MCD34B1/'
        landFile     = '/nobackup/'+inst+'/LevelB/invariant/'+inst.lower()+'-g5nr.lb2.asm_Nx.nc4'
    elif os.path.exists(nccs): 
        geoFile      = nccs  +inst.lower()+'.lg1.invariant.nc4'
        scanTimeFile = nccs + inst.lower()+'.lg1.scantime.nc4'
    else:
        geoFile      = inst.lower()+'.lg1.invariant.nc4'
        scanTimeFile = inst.lower()+'.lg1.scantime.nc4'	
        datadir      = '/nobackup/3/pcastell/MODIS/MCD43B1/'


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
        raise ValueError, 'invalid extension <%s>'%ext
    options.zlib = not options.nozip

    # Create (x,y,t) coordinates
    # --------------------------
    ncGeo, clon, clat, MISSING = getCoords(options.geoFile,options.verbose)
    nNS, nEW = clon.shape            

    
    for t in t1:
      path = datadir + t.strftime('%Y%j') + '/'
      assert os.path.exists(path), path + ' DOES NOT EXIST'
      mcdData = mcd43.McD43(path,clon[~clon.mask],clat[~clat.mask],Verb=0)

      # Gap Fill Data
      gap_fill(mcdData,clon,clat,options)

      writenc(mcdData,ncGeo,clon,clat,options)