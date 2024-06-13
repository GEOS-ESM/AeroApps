#!/usr/bin/env python3

"""
    Utility to sample G5NR model output in a SBG imager swath

    Adapted from tile_sampler.py
    Adapted from ACCP/lidar_sampler.py
    Trying to modernize to use xarray

    P. Castellanos Jun 2024
    P. Castellanos Feb 2024

"""

import os
import sys

import argparse
from   dateutil.parser import parse         as isoparser
from   datetime        import datetime, timedelta
from   pyobs.sgp4      import getTrack as getTrackTLE
import numpy           as np
from   MAPL            import eta
from   MAPL.constants  import *
from   pyobs.nc4ctl    import NC4ctl  
import xarray          as xr
from   netCDF4 import Dataset

class NC4ctl_(NC4ctl):
    interpXY = NC4ctl.interpXY_LatLon # select this as the default XY interpolation
#---

class TleVar(object):
    """
    Generic container for Variables
    """
    def __init__(self,name):
        self.name = name

#---
def Open(filename):
    """
    Uses GFIO or GFIOctl to open either a NetCDF-4 or a control file.
    Very heuristic.
    """
    from gfio import GFIOurl, GFIOctl
    name, ext = os.path.splitext(filename)

    # Open the GFIO dataset
    # ---------------------
    f = None
    if 'HTTP://' == name[:7].upper():
        f = GFIOurl(filename)
        f.lower = True # force variable names to be lower case when sampling.
    elif ext.upper() in ('.NC4','.NC','.HDF','.H5'):
        f = GFIOurl(filename)
        f.lower = False
    else:
        f = GFIOctl(filename)
        f.lower = False

    # Create variable dictionary
    # --------------------------
    Vars = dict()
    if len(f.vtitle)<len(f.vname):
        f.vtitle = f.vname[:]      # In case vtitle is not filled (hack)
    for i in range(len(f.vname)):
        if f.lower:
            v = f.vname[i].upper()
        else:
            v = f.vname[i]
        var = TleVar(v)
        var.title = f.vtitle[i]
        var.km = f.kmvar[i]
        if var.km>0:
            var.levunits = f.levunits[:]
            var.levs = f.levs[:]
        try:
            var.units = f.vunits[i]
        except:
            var.units = 'unknown'  # opendap currently lacks units
        Vars[v] = var

    f.Vars = Vars

    return f

def getVars(rcFile):
    """
    Parse reource file, create variable dictionary with relevant
    metadata.
    """
    from MAPL.config    import Config

    cf = Config(rcFile)
    Vars = dict()
    AllVars = dict()
    levUnits = 'none'
    levs = []
    for V in list(cf.keys()):
        path = cf(V)
        f = Open(path)
        varList = []
        if '*' in V:
            VARS = list(f.Vars.keys())
        else:
            VARS = V.split(',')
        for v in VARS:
            v = v.strip()
            var = f.Vars[v]
            if AllVars.__contains__(v): # unique variable names for output
                print(" >< Skipping duplicate variable <%s> in %s"%(v,path))
                continue
            elif v.upper() == "TAITIME":
                continue # annoying HDFEOS crap
            else:
                AllVars[v] = True
            if var.km>0:
                levUnits = var.levunits
                levs = var.levs
            varList += [var,]
        Vars[path] = varList

    return (Vars, levs, levUnits)


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


def writeNC ( args, Vars, levs, levUnits, 
              title='GEOS-5 Imager Sampler',
              doAkBk=False, zlib=True):
    """
    Write a NetCDF file with sampled GEOS-5 variables on a satellite swath
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

    km = len(levs)
    
    # Open NC file
    # ------------
    nc = Dataset(args.outFile,'w',format=args.format)

    # Set global attributes
    # ---------------------
    nc.title = title
    nc.institution = 'NASA/Goddard Space Flight Center'
    nc.source = 'Global Model and Assimilation Office'
    nc.history = 'Created from GEOS-5 standard collections by imager_sampler.py'
    nc.references = 'n/a'
    nc.comment = 'This file contains GEOS-5 related parameters along a satellite swath'
    nc.contact = 'Patricia Castellanos <patricia.castellanos@nasa.gov>'
    nc.Conventions = 'CF'
    nc.swathFile = args.swathFile
 
    # Create dimensions
    # -----------------
    nt = nc.createDimension('time',ds.sizes['time']) 
    ls = nc.createDimension('ls',int(str(ds.isotime.dtype)[-2:]))
    nalong = nc.createDimension('nalong',ds.sizes['nalong'])
    ncross = nc.createDimension('ncross',ds.sizes['ncross'])
    if km>0:
        nz = nc.createDimension('lev',km)
        if doAkBk:
            ne = nc.createDimension('ne',km+1)


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

    if km > 0: # pressure level not supported yet
        lev = nc.createVariable('lev','f4',('lev',),zlib=zlib)
        lev.long_name = 'Vertical Level'
        lev.units = levUnits.strip()
        lev.positive = 'down'
        lev.axis = 'z'
        lev[:] = levs[:]

        if doAkBk:
            ae, be = eta.getEdge(km) # Coefficients for Hybrid coordinates
            ak = nc.createVariable('ak','f4',('ne',),zlib=zlib)
            ak.long_name = 'Eta coordinate coefficient ak (p = ak + bk * ps)'
            ak.units = 'Pa'
            ak = ae[:]
            bk = nc.createVariable('bk','f4',('ne',),zlib=zlib)
            bk.long_name = 'Eta coordinate coefficient bk (p = ak + bk * ps)'
            bk.units = '1'
            bk = be[:]
    
    
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

      
    # Loop over datasets, sample and write each variable
    # --------------------------------------------------
    ntime,nalong,ncross  = ds.longitude.shape
    for path in Vars:
        if args.verbose:
            print(" <> opening "+path)
        g = Open(path) 
        g.nc4 = NC4ctl_(path)
        for var in Vars[path]:
            if var.km == 0:
                dim = ('time','nalong','ncross')
            else:
                dim = ('time','lev','nalong','ncross')
            this = nc.createVariable(var.name,'f4',dim,zlib=zlib)
            this.standard_name = var.title
            this.long_name = var.title.replace('_',' ')
            this.missing_value = np.float32(MAPL_UNDEF)
            this.units = var.units
            if g.lower:
                name = var.name.lower() # GDS always uses lower case
            else:
                name = var.name
            if args.verbose:
                print(" [] %s interpolating <%s>"%\
                        (args.algorithm.capitalize(),name.upper()))   

            # Use NC4ctl for linear interpolation
            # -----------------------------------
            for a in range(nalong):
                for c in range(ncross):
                    if g.tbeg == g.tend:
                        # this is a constant file
                        tconst = np.repeat(g.tbeg,ntime)
                        Z = g.nc4.sample(name,ds.isel(ncross=c,nalong=a).longitude.values,ds.isel(ncross=c,nalong=a).latitude.values,tconst,
                            Transpose=True,squeeze=True)
                    else:
                        Z = g.nc4.sample(name,ds.isel(ncross=c,nalong=a).longitude.values,ds.isel(ncross=c,nalong=a).latitude.values,tyme,
                             Transpose=True,squeeze=True)

                    Z[abs(Z)>MAPL_UNDEF/1000.] = MAPL_UNDEF # detect undef contaminated interp

                    if var.km > 0:
                        this[:,:,a,c] = Z
                    else:
                        this[:,a,c] = Z
            
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
    algo = 'linear'
    

#   Parse command line options
#   --------------------------
    parser = argparse.ArgumentParser()
    parser.add_argument("swathFile",
                        help="SBG swath file with lat/lons")

    parser.add_argument("rcFile",
                        help="model collection rcFile")

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
  
    # Get Variables and Metadata
    # --------------------------
    Vars, levs, levUnits = getVars(args.rcFile)

    # Write output file
    # -----------------
    writeNC(args,Vars,levs,levUnits)
    
