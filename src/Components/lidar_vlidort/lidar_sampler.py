#!/usr/bin/env python

"""
    Utility to generate a orbital trajectory and sample model
    output along it.

    Adapted from trj_sampler, but using NC4ctl to do interpolation
    as this is faster.

    P. Castellanos May 2017

"""

import os
import sys
if os.path.exists('/discover/nobackup'):
    sys.path.append(os.environ["HOME"]+'/workspace/GAAS/src/GMAO_Shared/GMAO_pyobs')
else:
    sys.path.append(os.environ["NOBACKUP"]+'/workspace/GAAS/src/GMAO_Shared/GMAO_pyobs')

from   trj_sampler     import getVars, getTrackICT, getTrackCSV, getTrackNPZ, getTrackHSRL, writeXLS
from   optparse        import OptionParser
from   dateutil.parser import parse         as isoparser
from   pyobs.sgp4      import getTrack as getTrackTLE
from   numpy           import zeros, arange, array
from   MAPL            import eta
from   MAPL.constants  import *
from   pyobs.nc4ctl    import NC4ctl  

class NC4ctl_(NC4ctl):
    interpXY = NC4ctl.interpXY_LatLon # select this as the default XY interpolation
#---
def writeNC ( lons, lats, tyme, Vars, levs, levUnits, trjFile, options,
              title='GEOS-5 Trajectory Sampler',
              doAkBk=False, zlib=False):
    """
    Write a NetCDF file with sampled GEOS-5 variables along the satellite track
    described by (lon,lat,tyme).
    """
    from netCDF4 import Dataset

    km = len(levs)
    
    # Open NC file
    # ------------
    nc = Dataset(options.outFile,'w',format=options.format)

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
    nc.trjFile = trjFile
 
    # Create dimensions
    # -----------------
    nt = nc.createDimension('time',len(tyme))
    ls = nc.createDimension('ls',19)
    if km>0:
        nz = nc.createDimension('lev',km)
        if doAkBk:
            ne = nc.createDimension('ne',km+1)
    x = nc.createDimension('x',1)
    y = nc.createDimension('y',1)

    # Coordinate variables
    # --------------------
    time = nc.createVariable('time','i4',('time',),zlib=zlib)
    time.long_name = 'Time'
    t0 = tyme[0]
    time.units = 'seconds since %s'%t0.isoformat(' ')
    time[:] = array([(t-t0).total_seconds() for t in tyme])
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
    
    # Add fake dimensions for GrADS compatibility
    # -------------------------------------------
    x = nc.createVariable('x','f4',('x',),zlib=zlib)
    x.long_name = 'Fake Longitude for GrADS Compatibility'
    x.units = 'degrees_east'
    x[:] = zeros(1)
    y = nc.createVariable('y','f4',('y',),zlib=zlib)
    y.long_name = 'Fake Latitude for GrADS Compatibility'
    y.units = 'degrees_north'
    y[:] = zeros(1)
    
    # Trajectory coordinates
    # ----------------------
    lon = nc.createVariable('trjLon','f4',('time',),zlib=zlib)
    lon.long_name = 'Trajectory Longitude'
    lon.units = 'degrees_east'
    lon[:] = lons[:]
    lat = nc.createVariable('trjLat','f4',('time',),zlib=zlib)
    lat.long_name = 'Trajectory Latitude'
    lat.units = 'degrees_north'
    lat[:] = lats[:]
    
    # Time in ISO format if so desired
    # ---------------------------------
    if options.isoTime:
        isotime = nc.createVariable('isotime','S1',('time','ls'),zlib=zlib)
        isotime.long_name = 'Time (ISO Format)'
        isotmp = zeros((len(lons),19),dtype='S1')
        for i in range(len(lons)):
            isotmp[i][:] = list(tyme[i].isoformat())
        isotime[:] = isotmp[:]
      
    # Loop over datasets, sample and write each variable
    # --------------------------------------------------
    for path in Vars:
        if options.verbose:
            print " <> opening "+path
        g = Open(path) 
        g.nc4 = NC4ctl_(path)
        for var in Vars[path]:
            if var.km == 0:
                dim = ('time',)
            else:
                dim = ('time','lev')
            this = nc.createVariable(var.name,'f4',dim,zlib=zlib)
            this.standard_name = var.title
            this.long_name = var.title.replace('_',' ')
            this.missing_value = MAPL_UNDEF
            this.units = var.units
            if g.lower:
                name = var.name.lower() # GDS always uses lower case
            else:
                name = var.name
            if options.verbose:
                print " [] %s interpolating <%s>"%\
                        (options.algo.capitalize(),name.upper())   

            # Use NC4ctl for linear interpolation
            # -----------------------------------
            Z = g.nc4.sample(name,lons,lats,tyme,
                             Transpose=True,squeeze=True)

            # Z = g.sample(name,lons,lats,tyme,algorithm=options.algo,
            #              Transpose=True,squeeze=True)
            Z[abs(Z)>MAPL_UNDEF/1000.] = MAPL_UNDEF # detect undef contaminated interp
            this[:] = Z
            
    # Close the file
    # --------------
    nc.close()

    if options.verbose:
        print " <> wrote %s file %s"%(options.format,options.outFile)
    
#---

#------------------------------------ M A I N ------------------------------------

if __name__ == "__main__":
    
    format = 'NETCDF3_CLASSIC'
    rcFile  = 'trj_sampler.rc'
    outFile = 'trj_sampler.nc'
    dt_secs = 60
    algo = 'linear'

#   Parse command line options
#   --------------------------
    parser = OptionParser(usage="Usage: %prog [OPTIONS] tleFile|ictFile|csvFile|npzFile [iso_t1 iso_t2]",
                          version='1.0.1' )

    parser.add_option("-a", "--algorithm", dest="algo", default=algo,
              help="Interpolation algorithm, one of linear, nearest (default=%s)"\
                          %algo)

    parser.add_option("-r", "--rcFile", dest="rcFile", default=rcFile,
              help="Resource file defining parameters to sample (default=%s)"\
                          %rcFile )

    parser.add_option("-o", "--output", dest="outFile", default=outFile,
              help="Output NetCDF file (default=%s)"\
                          %outFile )

    parser.add_option("-f", "--format", dest="format", default=format,
              help="Output file format: one of NETCDF4, NETCDF4_CLASSIC, NETCDF3_CLASSIC, NETCDF3_64BIT or EXCEL (default=%s)"%format )

    parser.add_option("-t", "--trajectory", dest="traj", default=None,
              help="Trajectory file format: one of tle, ict, csv, npz, hsrl (default=trjFile extension)" )

    parser.add_option("-d", "--dt_secs", dest="dt_secs", default=dt_secs,
              type='int',
              help="Timesetp in seconds for TLE sampling (default=%s)"%dt_secs )

    parser.add_option("-I", "--isoTime",
                      action="store_true", dest="isoTime",
                      help="Include ISO format time in output file.")

    parser.add_option("-v", "--verbose",
                      action="store_true", dest="verbose",
                      help="Verbose mode.")

    (options, args) = parser.parse_args()
    
    if len(args) == 3:
        trjFile, iso_t1, iso_t2 = args
        t1, t2 = (isoparser(iso_t1), isoparser(iso_t2))
    elif len(args) == 1:
        trjFile = args[0]
        t1, t2 = None, None
    else:
        parser.error("must have 1 or 3 arguments: tleFile|ictFile [iso_t1 iso_t2]")

    if options.traj is None:
        name, ext = os.path.splitext(trjFile)
        options.traj = ext[1:]
    options.traj = options.traj.upper()
        
    # Create consistent file name extension
    # -------------------------------------
    name, ext = os.path.splitext(options.outFile)
    if ext.upper() == '.XLS':
        options.format = 'EXCEL'
    if 'NETCDF4' in options.format:
        options.outFile = name + '.nc4'
    elif 'NETCDF3' in options.format:
        options.outFile = name + '.nc'
    elif 'EXCEL' in options.format:
        options.outFile = name + '.xls'
    else:
        raise ValueError, 'invalid extension <%s>'%ext
  
    # Get Variables and Metadata
    # --------------------------
    Vars, levs, levUnits = getVars(options.rcFile)

    # Create trajectory
    # -----------------
    if options.traj == 'TLE':
        if t1 is None:
            raise ValueError, 'time range (t1,t2) must be specified when doing TLE sampling.'
        lon, lat, tyme = getTrackTLE(trjFile, t1, t2, options.dt_secs)
    elif options.traj == 'ICT':
        lon, lat, tyme = getTrackICT(trjFile,options.dt_secs)
    elif options.traj == 'CSV':
        lon, lat, tyme = getTrackCSV(trjFile)
    elif options.traj == 'NPZ':
        lon, lat, tyme = getTrackNPZ(trjFile)
    elif options.traj == 'HSRL' or options.traj == 'H5':
        lon, lat, tyme = getTrackHSRL(trjFile)
    else:
        raise ValueError, 'cannot handle trajectory file format <%s>'%options.traj
    # Make sure longitudes in [-180,180]
    # ----------------------------------
    if lon.max()>180.:
        lon[lon>180] = lon[lon>180] - 360.
    
    # Write output file
    # -----------------
    if options.format == 'EXCEL':
        writeXLS(lon,lat,tyme,Vars,trjFile,options)
    else:
        writeNC(lon,lat,tyme,Vars,levs,levUnits,trjFile,options)
    
