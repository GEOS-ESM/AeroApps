#!/usr/bin/env python3

"""
    Utility to sample G5NR model output for 1 pixel of an SBG tile 

    For doing quick and dirty single pixel testing

    Adapted from ACCP/lidar_sampler.py
    Trying to modernize to use xarray

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


def writeNC ( args, tyme, Vars, levs, levUnits, 
              title='GEOS-5 Trajectory Sampler',
              doAkBk=False, zlib=True):
    """
    Write a NetCDF file with sampled GEOS-5 variables on a satellite tile
    described by (lon,lat,tyme).
    """

    # read lons/lats
    # ----------------
    ds = xr.open_dataset(args.tileFile,group='ancillary')

    # Make sure longitudes in [-180,180]
    # ----------------------------------
    if ds.lon.max()>180.:
        ds.lon[ds.lon>180] = ds.lon[ds.lon>180] - 360.


    from netCDF4 import Dataset

    km = len(levs)
    
    # Open NC file
    # ------------
    nc = Dataset(args.outFile,'w',format=args.format)

    # Set global attributes
    # ---------------------
    nc.title = title
    nc.institution = 'NASA/Goddard Space Flight Center'
    nc.source = 'Global Model and Assimilation Office'
    nc.history = 'Created from GEOS-5 standard collections by trj_sampler.py'
    nc.references = 'n/a'
    nc.comment = 'This file contains GEOS-5 related parameters along a satellite or aircraft track.'
    nc.contact = 'Arlindo da Silva <arlindo.dasilva@nasa.gov>'
    nc.Conventions = 'CF'
    nc.tileFile = args.tileFile
 
    # Create dimensions
    # -----------------
    nt = nc.createDimension('time',None) 
    ls = nc.createDimension('ls',19)
    x = nc.createDimension('lon',len(ds.lon))
    y = nc.createDimension('lat',len(ds.lat))
    if km>0:
        nz = nc.createDimension('lev',km)
        if doAkBk:
            ne = nc.createDimension('ne',km+1)


    # Coordinate variables
    # --------------------
    time = nc.createVariable('time','i4',('time',),zlib=zlib)
    time.long_name = 'Time'
    t0 = tyme[0]
    time.units = 'seconds since %s'%t0.isoformat(' ')
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
    lon = nc.createVariable('trjLon','f4',('lon',),zlib=zlib)
    lon.long_name = 'Tile Longitude'
    lon.units = 'degrees_east'
    lon[:] = ds.lon
    lat = nc.createVariable('trjLat','f4',('lat',),zlib=zlib)
    lat.long_name = 'Tile Latitude'
    lat.units = 'degrees_north'
    lat[:] = ds.lat
    
    # Time in ISO format if so desired
    # ---------------------------------
    if args.isoTime:
        isotime = nc.createVariable('isotime','S1',('time','ls'),zlib=zlib)
        isotime.long_name = 'Time (ISO Format)'
        isotmp = np.zeros((len(tyme),19),dtype='S1')
        for i in range(len(tyme)):
            isotmp[i][:] = list(tyme[i].isoformat())
        isotime[:] = isotmp[:]
      
    # Loop over datasets, sample and write each variable
    # --------------------------------------------------
    glon,glat = np.meshgrid(ds.lon,ds.lat)
    nlon,nlat = len(ds.lon),len(ds.lat)
    for path in Vars:
        if args.verbose:
            print(" <> opening "+path)
        g = Open(path) 
        g.nc4 = NC4ctl_(path)
        for var in Vars[path]:
            if var.km == 0:
                dim = ('time','lat','lon')
            else:
                dim = ('time','lat','lon','lev')
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
                print(" [] Linear interpolating <%s>"%\
                        (name.upper()))   

            # Use NC4ctl for linear interpolation
            # -----------------------------------
            lons = np.array([glon[args.row,args.row]])
            lats = np.array([glat[args.row,args.col]])
            z = g.nc4.sample(name,lons,lats,tyme,
                             Transpose=True,squeeze=True)

            z[abs(z)>MAPL_UNDEF/1000.] = MAPL_UNDEF # detect undef contaminated interp

            if var.km == 0:
                Z = np.ones([1,nlat,nlon])*MAPL_UNDEF
                Z[0,args.row,args.col] = z
            else:
                Z = np.ones([1,nlat,nlon,var.km])
                Z[0,args.row,args.col,:] = z
            
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
    
    format = 'NETCDF4'
    outFile = 'tile_sampler.nc'
    

#   Parse command line options
#   --------------------------
    parser = argparse.ArgumentParser()
    parser.add_argument("tileFile",
                        help="SBG tile file with lat/lons")


    parser.add_argument("rcFile",
                        help="model collection rcFile")

    parser.add_argument("row",type=int,
                        help="row index")

    parser.add_argument("col",type=int,
                        help="col index")

    parser.add_argument("-o", "--outFile", default=outFile,
              help="Output NetCDF file (default=%s)"\
                          %outFile )

    parser.add_argument("-f", "--format", default=format,
              help="Output file format: one of NETCDF4, NETCDF4_CLASSIC, NETCDF3_CLASSIC, NETCDF3_64BIT(default=%s)"%format )

    parser.add_argument("--const", action="store_true",
              help="This is a constant file" )

    parser.add_argument("-I", "--isoTime", action="store_true", 
                      help="Include ISO format time in output file.")

    parser.add_argument("-v", "--verbose",
                      action="store_true", dest="verbose",
                      help="Verbose mode.")

    args = parser.parse_args()

    if args.const:
        tyme = np.array([datetime(2005,5,16,0)])
    else:
        # get julian day from file name
        jday = args.tileFile.split('_')[-2][4:]

        # need lon to get the UTC time
        # use the center lon
        # assume satellite passes at 1300 local solar time
        with  xr.open_dataset(args.tileFile, group="ancillary") as ds:
            lon0 = float(ds.lon.mean()) 

        dmin = int(60.*lon0/15.)
        
        # get date
        tyme  = datetime.strptime('2006-{}T13'.format(jday),'%Y-%jT%H')
        tyme  -= timedelta(minutes=dmin)
        tyme  = np.array([tyme])

    
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
    writeNC(args,tyme,Vars,levs,levUnits)
    
