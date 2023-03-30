#!/usr/bin/env python3

"""
    Utility to generate a orbital trajectory and sample model
    output along it. 

"""

import os

from numpy import zeros, arange, array

from optparse        import OptionParser
from datetime        import datetime, timedelta
from dateutil.parser import parse         as isoparser

from MAPL           import Config, eta
from MAPL.constants import *
from sgp4           import getTrack

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
    for i in range(len(f.vname)):
        v = f.vname[i].upper()
        var = TleVar(v)
        var.title = f.vtitle[i]
        var.km = f.kmvar[i]
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

    cf = Config(rcFile)
    Vars = dict()
    km = 0
    for V in list(cf.keys()):
        path = cf(V)
        f = Open(path)
        varList = []
        for v in V.split(','):
            var = f.Vars[v]
            km = max(km,var.km)
            varList += [var,]
        Vars[path] = varList
        
    return (Vars, km)

#---
def writeNC ( lons, lats, tyme, Vars, km, filename, tle,
              format='NETCDF3', zlib=False, isoTime=True,
              title='GEOS-5 TLE Sampler', Verbose=False):
    """
    Write a NetCDF file with sampled GEOS-5 variables along the satellite track
    described by (lon,lat,tyme).
    """
    from netCDF4 import Dataset

    # Open NC file
    # ------------
    nc = Dataset(filename,'w',format=format)

    # Set global attributes
    # ---------------------
    nc.title = title
    nc.institution = 'NASA'
    nc.source = 'Global Model and Assimilation Office'
    nc.history = 'Created from GEOS-5 standard collections by tle_sampler.py'
    nc.references = 'n/a'
    nc.comment = 'This file contains GEOS-5 cloud related parameters along a satellite track.'
    nc.contact = 'Arlindo da Silva <arlindo.dasilva@nasa.gov>'
    nc.Conventions = 'CF'
    nc.TLE = tle
 
    # Create dimensions
    # -----------------
    nz = nc.createDimension('nz',km)
    ne = nc.createDimension('ne',km+1)
    nt = nc.createDimension('nt',len(tyme))
    ls = nc.createDimension('ls',19)

    # Coordinate variables
    # --------------------
    lon = nc.createVariable('lon','f4',('nt',),zlib=zlib)
    lon.long_name = 'Longitude'
    lon.units = 'degrees_east'
    lon[:] = lons[:]
    
    lat = nc.createVariable('lat','f4',('nt',),zlib=zlib)
    lat.long_name = 'Latitude'
    lat.units = 'degrees_north'
    lat[:] = lats[:]

    lev = nc.createVariable('lev','f4',('nz',),zlib=zlib)
    lev.long_name = 'Model Level for Hybrid Sigma-pressure Coordinates'
    lev.units = '1'
    lev.positive = 'down'
    lev.axis = 'z'
    lev[:] = arange(1,km+1)

    ae, be = eta.getEdge(km) # Coefficients for Hybrid coordinates

    ak = nc.createVariable('ak','f4',('ne',),zlib=zlib)
    ak.long_name = 'Eta coordinate coefficient ak (p = ak + bk * ps)'
    ak.units = 'Pa'
    ak = ae[:]
    
    bk = nc.createVariable('bk','f4',('ne',),zlib=zlib)
    bk.long_name = 'Eta coordinate coefficient bk (p = ak + bk * ps)'
    bk.units = '1'
    bk = be[:]
    
    time = nc.createVariable('time','i4',('nt',),zlib=zlib)
    time.long_name = 'Time'
    t0 = tyme[0]
    time.units = 'seconds since %s'%t0.isoformat(' ')
    time[:] = array([(t-t0).total_seconds() for t in tyme])

    # Time in ISO format if so desired
    # ---------------------------------
    if isoTime:
        isotime = nc.createVariable('isotime','S1',('nt','ls'),zlib=zlib)
        isotime.long_name = 'Time (ISO Format)'
        isotmp = zeros((len(lons),19),dtype='S1')
        for i in range(len(lons)):
            isotmp[i][:] = list(tyme[i].isoformat())
        isotime[:] = isotmp[:]
      
    # Loop over datasets, sample and write each variable
    # --------------------------------------------------
    for path in Vars:
        if Verbose:
            print(" <> opening "+path)
        g = Open(path) 
        for var in Vars[path]:
            if var.km == 0:
                dim = ('nt',)
            else:
                dim = ('nt','nz')
            this = nc.createVariable(var.name,'f4',dim,zlib=zlib)
            this.standard_name = var.title
            this.long_name = var.title.replace('_',' ')
            this.missing_value = MAPL_UNDEF
            this.units = var.units
            if g.lower:
                name = lower(var.name) # GDS always uses lower case
            else:
                name = var.name
                    
            this_ = g.sample(name,lons,lats,tyme,Verbose=Verbose)
            this[:] = this_.T # transpose it so that shape is (nt,nz)
            
    # Close the file
    # --------------
    nc.close()

    if Verbose:
        print(" <> wrote %s file %s"%(format,filename))
    
#------------------------------------ M A I N ------------------------------------

if __name__ == "__main__":
    
    format = 'NETCDF3_CLASSIC'
    rcFile  = 'tle_sampler.rc'
    outFile = 'tle_sampler.nc'
    dt_secs = 60

#   Parse command line options
#   --------------------------
    parser = OptionParser(usage="Usage: %prog [OPTIONS] tleFile iso_t1 iso_t2",
                          version='1.0.0' )

    parser.add_option("-r", "--rcFile", dest="rcFile", default=rcFile,
              help="Resource file defining parameters to sample (default=%s)"\
                          %rcFile )

    parser.add_option("-o", "--output", dest="outFile", default=outFile,
              help="Output NetCDF file (default=%s)"\
                          %outFile )

    parser.add_option("-f", "--format", dest="format", default=format,
              help="Output file format: one of NETCDF4, NETCDF4_CLASSIC, NETCDF3_CLASSIC or NETCDF3_64BIT (default=%s)"%format )

    parser.add_option("-d", "--dt_secs", dest="dt_secs", default=dt_secs,
              type='int',
              help="Timesetp in seconds (default=%s)"%dt_secs )

    parser.add_option("-I", "--isoTime",
                      action="store_true", dest="isoTime",
                      help="Include time in ISO format as well.")

    parser.add_option("-v", "--verbose",
                      action="store_true", dest="verbose",
                      help="Verbose mode.")


    (options, args) = parser.parse_args()
    
    if len(args) == 3:
        tleFile, iso_t1, iso_t2 = args
        t1, t2 = (isoparser(iso_t1), isoparser(iso_t2))
    else:
        parser.error("must have 3 arguments: tleFile iso_t1 iso_t2")

    # Create consistent file name extension
    # -------------------------------------
    name, ext = os.path.splitext(options.outFile)
    if 'NETCDF4' in options.format:
        options.outFile = name + '.nc4'
    else:
        options.outFile = name + '.nc'

    # Get Variables and Metadata
    # --------------------------
    Vars, km = getVars(options.rcFile)

    # Create trajectory
    # -----------------
    lon, lat, tyme = getTrack(tleFile, t1, t2, options.dt_secs)

    # Write output file
    # -----------------
    writeNC(lon,lat,tyme,Vars,km,options.outFile,tleFile,
            format=options.format, isoTime=options.isoTime,
            Verbose=options.verbose)

    
