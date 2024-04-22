#!/usr/bin/env python3

"""
    Utility to sample G5NR model output in a SBG tile

    Adapted from ACCP/lidar_sampler.py
    Trying to modernize to use xarray

    P. Castellanos Feb 2024

"""

import os
import sys

import argparse
from   dateutil.parser import parse         as isoparser
from   datetime        import datetime, timedelta
import numpy           as np
from   MAPL            import eta
from   MAPL.config     import Config
from   MAPL.constants  import *
from   pyobs.nc4ctl    import NC4ctl  
import xarray          as xr
from   glob            import glob


def writeNC ( args, tyme, Vars, levs, levUnits, 
              title='GEOS-5 Trajectory Sampler',
              doAkBk=False, zlib=False):
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
        isotmp = zeros((len(tyme),19),dtype='S1')
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
                print(" [] %s interpolating <%s>"%\
                        (args.algo.capitalize(),name.upper()))   

            # Use NC4ctl for linear interpolation
            # -----------------------------------
            Z = g.nc4.sample(name,glon.ravel(),glat.ravel(),tyme.repeat(nlon*nlat),
                             Transpose=True,squeeze=True)

            Z[abs(Z)>MAPL_UNDEF/1000.] = MAPL_UNDEF # detect undef contaminated interp

            if var.km == 0:
                Z = Z.reshape([1,nlat,nlon])
            else:
                Z = Z.reshape([1,nlat,nlon,var.km])
            
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
    
    algo   = 'linear'
    dt_days = 1

#   Parse command line options
#   --------------------------
    parser = argparse.ArgumentParser()
    parser.add_argument("tile",
                        help="SBG tile")

    parser.add_argument("orbit",
                        help="orbit name")

    parser.add_argument("col",
                        help="collection name")

    parser.add_argument("swath_pcf",
                        help="swath filenames pcf input file")

    parser.add_argument("tile_pcf",
                        help="tile filenames pcf input file")

    parser.add_argument("start_isotime",
                        help="start isotime")

    parser.add_argument("end_isotime",
                        help="end isotime")

    parser.add_argument("-a", "--algorithm", default=algo,
              help="Interpolation algorithm, one of linear, nearest (default=%s)"\
                          %algo)

    parser.add_argument("-d", "--dt_days", default=dt_days,type=int,
              help="Timesetp in hours for swath files (default=%s)"%dt_days )

    parser.add_argument("-v", "--verbose",
                      action="store_true", dest="verbose",
                      help="Verbose mode.")

    args = parser.parse_args()

    cf_swath = Config(args.swath_pcf,delim=' = ')
    cf_tile  = Config(args.tile_pcf,delim=' = ')

    # get lat/lon for tile 
    nlon,nlat = int(cf_tile('nlon')), int(cf_tile('nlat'))
    H,V = int(args.tile[1:3]),int(args.tile[4:])

    lats = 60.0 - (V*6 + np.arange(nlat)*0.01) - 0.005
    lons = -180.0 + (H*6 + np.arange(nlon)*0.01) + 0.005

    # loop through start to end time
    sdate = isoparser(args.start_isotime)
    edate = isoparser(args.end_isotime)
    dt    = timedelta(days=args.dt_days)


    while sdate < edate:
      
        year = sdate.strftime('%Y')
        month = sdate.strftime('%m')
        day   = sdate.strftime('%d')
        nymd  = year+month+day

        # Get Filenames 
        # -------------
        inDir = cf_swath('inDir').replace('%ORBITNAME',args.orbit.upper()).replace('%year',year).replace('%month',month).replace('%day',day)
        for hh in range(24):
            swathFile = inDir + '/' + cf_swath('inFile').replace('%orbitname',args.orbit.lower()).replace('%col',args.col).replace('%nymd',nymd).replace('%hour',"{:02d}".format(hh))

        swath = xr.open_mfdataset(swathList)

        # Find the part of the swath that intersect with the tile
        # ---------------------------------------
        


#        # Write output file
#        # -----------------
#        writeNC(args,tyme,Vars,levs,levUnits)
    

        sdate += dt
