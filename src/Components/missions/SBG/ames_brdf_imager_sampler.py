#!/usr/bin/env python3

"""
    Utility to sample AMES brdf on a SBG imager swath

    Adapted from g5nr_imager_sampler.py
    Adapted from tile_sampler.py
    Adapted from ACCP/lidar_sampler.py
    Trying to modernize to use xarray

    P. Castellanos Jul 2024
    P. Castellanos Jun 2024
    P. Castellanos Feb 2024
"""

import os
import sys

import argparse
from   dateutil.parser import parse         as isoparser
from   datetime        import datetime, timedelta
import numpy           as np
from   MAPL.constants  import *
import xarray          as xr
from   netCDF4 import Dataset
from   ames_brdf import tiles, AMES_BRDF

SDS = ['Kg','Ki','Kv']

#---

def to_datetime(date):
    """
    Converts a numpy datetime64 object to a python datetime object 
    Input:
      date - a np.datetime64 object
    Output:
      DATE - a python datetime object
    """
    timestamp = ((date - np.datetime64('1970-01-01T00:00:00'))/ np.timedelta64(1, 's'))

    tyme = [datetime.utcfromtimestamp(tt) for tt in timestamp]
    return np.array(tyme)


def writeNC ( args, brdf, 
              title='AMES BRDF Imager Sampler',
              zlib=True):
    """
    Write a NetCDF file with sampled BRDF variables on a satellite swath
    described by (lon,lat,tyme).
    """

    # read lons/lats & tyme
    # ----------------
    ds = xr.open_dataset(args.swathFile)
    
    # convert numpy datetime to datetime
    tyme = to_datetime(ds.time.values)

    # get time units
    nc = Dataset(args.swathFile)
    dt_units = nc.variables['time'].units
    nc.close()

    # Make sure longitudes in [-180,180]
    # ----------------------------------
    lon = ds.longitude.values  # this is a pointer - changes to lon will occur in ds.longitude too, allows me to do numpy type indexing
    if lon.max()>180.:
        lon[lon>180] = lon[lon>180] - 360
    if lon.min()<-180.:
        lon[lon<-180] = lon[lon<-180] + 360

    # Open NC file
    # ------------
    nc = Dataset(args.outFile,'w',format=args.format)

    # Set global attributes
    # ---------------------
    nc.title = title
    nc.institution = 'NASA/Goddard Space Flight Center'
    nc.source = 'Global Model and Assimilation Office'
    nc.history = 'Created from GEOS-5 standard collections by ames_brdf_imager_sampler.py'
    nc.references = 'n/a'
    nc.comment = 'This file contains BRDF parameters along a satellite swath'
    nc.contact = 'Patricia Castellanos <patricia.castellanos@nasa.gov>'
    nc.Conventions = 'CF'
    nc.swathFile = args.swathFile
 
    # Create dimensions
    # -----------------
    nt = nc.createDimension('time',ds.sizes['time']) 
    ls = nc.createDimension('ls',int(str(ds.isotime.dtype)[-2:]))
    nalong = nc.createDimension('nalong',ds.sizes['nalong'])
    ncross = nc.createDimension('ncross',ds.sizes['ncross'])
    nwav   = nc.createDimension('nwav',brdf.nwav)

    # Coordinate variables
    # --------------------
    time = nc.createVariable('time','i4',('time',),zlib=zlib)
    time.long_name = 'Time'
    t0 = tyme[0]
    time.units = dt_units
    if 'microseconds' in dt_units:
        time[:] = np.array([(t-t0).total_seconds()*1e6 for t in tyme])
    else:
        time[:] = np.array([(t-t0).total_seconds() for t in tyme])

    wav = nc.createVariable('nwav','f4',('nwav',),zlib=zlib)
    wav.long_name = 'wavelength'
    wav.units = "micron (1e-6 meter)" 
    wav[:] = brdf.wavelength

    # Trajectory coordinates
    # ----------------------
    lon = nc.createVariable('trjLon','f4',('time',),zlib=zlib)
    lon.long_name = 'Trajectory Longitude'
    lon.units = 'degrees_east'
    lon[:] = ds.trjLon.values
    lat = nc.createVariable('trjLat','f4',('time',),zlib=zlib)
    lat.long_name = 'Trajectory Latitude'
    lat.units = 'degrees_north'
    lat[:] = ds.trjLat.values
    
    # Time in ISO format if so desired
    # ---------------------------------
    if args.isoTime:
        isotime = nc.createVariable('isotime','S1',('time','ls'),zlib=zlib)
        isotime.long_name = 'Time (ISO Format)'
        iso = [np.fromiter(t.decode(),(np.compat.unicode,1)) for t in ds.isotime.values]
        isotime[:] = iso

    # swath cooredinates
    lon = nc.createVariable('longitude','f4',('time','nalong','ncross'),zlib=zlib)
    lon.long_name = 'longitude of granule pixel centers'
    lon.units     = 'degrees_east'
    lon[:] = ds.longitude.values

    lat = nc.createVariable('latitude','f4',('time','nalong','ncross'),zlib=zlib)
    lat.long_name = 'latitude of granule pixel centers'
    lat.units     = 'degrees_north'
    lat[:] = ds.latitude.values


    # sample the brdf dataset to new lat/lon
    # --------------------------------------
    brdf.interp = brdf.ds.interp(lon=ds.longitude,
                                 lat=ds.latitude,
                                 method=args.algorithm)
      
    # Loop over variables and write to file
    # --------------------------------------------------
    ntime,nalong,ncross  = ds.longitude.shape
    for sds in SDS:
        if args.verbose:
            print(" [] %s writing <%s>"%(args.algorithm.capitalize(),sds))
        dim = ('nwav','time','nalong','ncross')
        this = nc.createVariable(sds,'f4',dim,zlib=zlib)
        this.standard_name = sds
        this.long_name = "RTLS BRDF Coefficient {}".format(sds)
        this.missing_value = np.float32(MAPL_UNDEF)
        this.units = "None"


        # Check for nan's - these are pixels outside of the brdf domain
        Z = brdf.interp[sds].to_numpy()
        I = np.isnan(Z)
        Z[I] = MAPL_UNDEF
        # Check the BRDF mask - this is where there is no data.  Usually ocean pixels 
        I = brdf.interp.mask == 0
        Z[:,I] = MAPL_UNDEF

        this[:] = Z
            
    # Close the file
    # --------------
    nc.close()
    ds.close()

    if args.verbose:
        print(" <> wrote %s file %s"%(args.format,args.outFile))
    
#---

#------------------------------------ M A I N ------------------------------------

if __name__ == "__main__":
    
    format = 'NETCDF4_CLASSIC'
    outFile = 'imager_sampler.nc4'
    dt_days = 8
    algo = 'nearest'
    

#   Parse command line options
#   --------------------------
    parser = argparse.ArgumentParser()
    parser.add_argument("swathFile",
                        help="SBG swath file with lat/lons")

    parser.add_argument("brdfPath",
                        help="path for the brdf dataset")

    parser.add_argument("-a", "--algorithm", default=algo,
              help="Interpolation algorithm, one of linear, nearest (default=%s)"\
                          %algo)

    parser.add_argument("-o", "--outFile", default=outFile,
              help="Output NetCDF file (default=%s)"\
                          %outFile )

    parser.add_argument("-f", "--format", default=format,
              help="Output file format: one of NETCDF4, NETCDF4_CLASSIC, NETCDF3_CLASSIC, NETCDF3_64BIT(default=%s)"%format )

    parser.add_argument("-I", "--isoTime", action="store_true", 
                      help="Include ISO format time in output file.")

    parser.add_argument("-v", "--verbose",
                      action="store_true", dest="verbose",
                      help="Verbose mode.")

    args = parser.parse_args()

    
    # Create consistent file name extension
    # -------------------------------------
    name, ext = os.path.splitext(args.outFile)
    if 'NETCDF4' in args.format:
        args.outFile = name + '.nc4'
    elif 'NETCDF3' in args.format:
        args.outFile = name + '.nc'
    else:
        raise ValueError('invalid extension <%s>'%ext)
  
    # figure out closest DOY to swath
    # ------------------------------
    ds = xr.open_dataset(args.swathFile)
    doy = int(to_datetime(np.array([ds.time.mean().values]))[0].strftime('%j'))
    doy = round(doy/8)*8 + 1
    doy = '{:03d}'.format(doy)

    # figure out the bounding box of the swath
    # ----------------------------------------
    box = [float(ds.longitude.min()),float(ds.longitude.max()),float(ds.latitude.min()),float(ds.latitude.max())]
    
    ds.close()

    # get tiles that overlap swath
    # ----------------------------
    grid = tiles(box)

    # if there are any tiles that overlap
    # read in brdf data
    # Write output file
    # ------------------
    if np.sum(grid) > 0:
        brdf = AMES_BRDF(args.brdfPath,doy,grid=grid,return_ds=True)
        writeNC(args,brdf)
    
