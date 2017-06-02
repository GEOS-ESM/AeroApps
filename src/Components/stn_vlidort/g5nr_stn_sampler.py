#!/usr/bin/env python

"""
    Utility to sample GEOS-5  files at fixed station locations.
    Adapted from stn_sampler to deal with G5NR file conventions

"""

import os
import sys
if os.path.exists('/discover/nobackup'):
    sys.path.append(os.environ["NOBACKUP"]+'/workspace/GAAS/src/GMAO_Shared/GMAO_pyobs')
else:
    sys.path.append(os.environ["HOME"]+'/workspace/GAAS/src/GMAO_Shared/GMAO_pyobs')
from trj_sampler  import getVars
from stn_sampler  import StnVar, Open, getStations
from numpy import zeros, ones, arange, array

from optparse        import OptionParser
from datetime        import datetime, timedelta
from dateutil.parser import parse         as isoparser
from csv             import DictReader

from MAPL           import Config
from MAPL.constants import *
from   pyobs.nc4ctl    import NC4ctl


class NC4ctl_(NC4ctl):
    interpXY = NC4ctl.interpXY_LatLon # select this as the default XY interpolation

def getTyme(Dt,t1,t2):
    t = t1
    tyme = [t,]
    while t < t2-Dt:
        t += Dt
        tyme += [t,]
    tyme = array(tyme)
    return tyme

#---
def stnSample(f,V,stnLon,stnLat,tyme,options,squeeze=True):
    """
    Sample file at station locations.
    """
    if options.verbose:
        print " [] Interpolating <%s>"%V.name
    ns, nt, nz = (len(stnLon), len(tyme), V.km)
    if f.lower:
        name = V.name.lower() # GDS always uses lower case
    else:
        name = V.name
    if nz>0:
        Z = zeros((ns,nt,nz))
    else:
        Z = zeros((ns,nt))
    n = 0
    for t in tyme:
        try:
            z = f.nc4.interpXY(name,stnLon,stnLat,t,algorithm=options.algo,
                         Transpose=True,squeeze=squeeze)            
            # z = f.interp(name,stnLon,stnLat,tyme=t,algorithm=options.algo,
            #              Transpose=True,squeeze=squeeze)
        except:
            print "    - Interpolation failed for <%s> on %s"%(V.name,str(t))
            if nz>0:
                z = MAPL_UNDEF * ones((ns,nz))
            else:
                z = MAPL_UNDEF * ones(ns)

        if nz>0:
            Z[:,n,:] = z
        else:
            Z[:,n] = z
        n += 1
        Z[abs(Z)>MAPL_UNDEF/1000.] = MAPL_UNDEF # detect undef contaminated interp
    return Z

#---
def writeNC ( stnName, stnLon, stnLat, tyme, Vars, levs, levUnits, options,
              zlib=False, title='GEOS-5 Station Sampler'):
    """
    Write a NetCDF file with sampled GEOS-5 variables at the station locations
    described by (lon,lat).
    """
    from netCDF4 import Dataset

    ns_, nt_, nz_ = ( len(stnLon), len(tyme), len(levs) )

    # Open NC file
    # ------------
    nc = Dataset(options.outFile,'w',format=options.format)

    # Set global attributes
    # ---------------------
    nc.Title = title
    nc.Institution = 'NASA/Goddard Space Flight Center'
    nc.Source = 'Global Model and Assimilation Office'
    nc.History = 'Created from GEOS-5 standard collections by stn_sampler.py'
    nc.References = 'n/a'
    nc.Comment = 'This file contains GEOS-5 parameters sampled at station locations.'
    nc.Contact = 'Arlindo da Silva <arlindo.dasilva@nasa.gov>'
    nc.Conventions = 'CF'

    # Create dimensions
    # -----------------
    ns = nc.createDimension('station', ns_ )
    if nz_>0:
        nz = nc.createDimension('lev', nz_ )
    nt = nc.createDimension('time', nt_ )
    ls = nc.createDimension('ls',19)
    x = nc.createDimension('x',1)
    y = nc.createDimension('y',1)
     
    # Station names
    # -------------
    stnName_ = nc.createVariable('stnName','S1',('station','ls',),zlib=zlib)
    stnName_.long_name = 'Station Names'
    stnName_.axis = 'e'
    stntmp = zeros((ns_,19),dtype='S1')
    for i in range(ns_):
        stntmp[i][:] = list('%-19s'%stnName[i])
    stnName_[:] = stntmp[:]

    # Coordinate variables
    # --------------------
    lon = nc.createVariable('stnLon','f4',('station',),zlib=zlib)
    lon.long_name = 'Longitude'
    lon.units = 'degrees_east'
    lon[:] = stnLon[:]
    
    lat = nc.createVariable('stnLat','f4',('station',),zlib=zlib)
    lat.long_name = 'Latitude'
    lat.units = 'degrees_north'
    lat[:] = stnLat[:]

    time = nc.createVariable('time','i4',('time',),zlib=zlib)
    time.long_name = 'Time'
    t0 = tyme[0]
    time.units = 'seconds since %s'%t0.isoformat(' ')
    time[:] = array([(t-t0).total_seconds() for t in tyme])    

    # Add fake dimensions for GrADS compatibility
    # ------------------------------------------
    x_ = nc.createVariable('x','f4',('x',),zlib=zlib)
    x_.long_name = 'Fake Longitude for GrADS Compatibility'
    x_.units = 'degrees_east'
    x_[:] = zeros(1)

    y_ = nc.createVariable('y','f4',('y',),zlib=zlib)
    y_.long_name = 'Fake Latitude for GrADS Compatibility'
    y_.units = 'degrees_north'
    y_[:] = zeros(1)

    en = nc.createVariable('station','i4',('station',),zlib=zlib)
    # e.long_name = 'Station Ensemble Dimension'
    # e.axis = 'e'
    # e.grads_dim = 'e'
    # e[:] = range(ns_)
    
    if nz_ > 0:
        lev = nc.createVariable('lev','f4',('lev',),zlib=zlib)
        lev.long_name = 'Vertical Level'
        lev.units = levUnits.rstrip()
        lev.positive = 'down'
        lev.axis = 'z'
        lev[:] = levs[:]


    # Time in ISO format if so desired
    # ---------------------------------
    if options.isoTime:
        isotime = nc.createVariable('isotime','S1',('time','ls',),zlib=zlib)
        isotime.long_name = 'Time (ISO Format)'
        isotmp = zeros((nt_,19),dtype='S1')
        for i in range(nt_):
            isotmp[i][:] = list(tyme[i].isoformat())
        isotime[:] = isotmp[:]
      
    # # Loop over variables on GFIO file.
    # # --------------------------------
    # for path in Vars:
    #     if options.verbose:
    #         print " <> opening "+path
    #     f = Open(path) 
    #     f.nc4 = NC4ctl_(path)
    #     for var in Vars[path]:
    #         if var.km == 0:
    #             dim = ('station','time',)
    #             shp = ( ns_, nt_)
    #         else:
    #             dim = ('station','time','lev',)
    #             shp = ( ns_, nt_, nz_)
    #         this = nc.createVariable(var.name,'f4',dim,zlib=zlib)
    #         this.standard_name = var.title
    #         this.long_name = var.title.replace('_',' ')
    #         this.missing_value = MAPL_UNDEF
    #         this.units = var.units
    #         if options.dryrun:
    #             this_ = zeros(shp)
    #             if options.verbose:
    #                 print "[] Zero-filling <%s>"%var.name
    #         else:
    #             this_ = stnSample(f,var,stnLon,stnLat,tyme,options)
    #         this[:] = this_[:]
            
    # Close the file
    # --------------
    nc.close()

    if options.verbose:
        print " <> wrote %s file %s"%(options.format,options.outFile)
    
#------------------------------------ M A I N ------------------------------------

if __name__ == "__main__":
    
    format = 'NETCDF3_CLASSIC'
    outFile = 'stn_sampler.nc'
    algo = 'linear'
    dt_secs = 60*30  #30 minutes
    
#   Parse command line options
#   --------------------------
    parser = OptionParser(usage="Usage: %prog [OPTIONS] stnFile g5File [iso_t1 iso_t2]",
                          version='1.0.0' )

    parser.add_option("-o", "--output", dest="outFile", default=outFile,
              help="Output NetCDF file (default=%s)"\
                          %outFile )

    parser.add_option("-a", "--algorithm", dest="algo", default=algo,
              help="Interpolation algorithm, one of linear, nearest (default=%s)"\
                          %algo)

    parser.add_option("-V", "--vars", dest="Vars", default=None,
              help="Variables to sample (default=All)")
    
    parser.add_option("-f", "--format", dest="format", default=format,
              help="Output file format: one of NETCDF4, NETCDF4_CLASSIC, NETCDF3_CLASSIC, NETCDF3_64BIT (default=%s)"%format )

    parser.add_option("-d", "--dt_secs", dest="dt_secs", default=dt_secs,
              type='int',
              help="Timesetp in seconds for sampling (default=%s)"%dt_secs )    

    parser.add_option("-I", "--isoTime",
                      action="store_true", dest="isoTime",
                      help="Include time in ISO format as well.")

    parser.add_option("-v", "--verbose",
                      action="store_true", dest="verbose",
                      help="Verbose mode.")

    parser.add_option("-n", "--dryrun",
                      action="store_true", dest="dryrun",
                      help="Dry-run mode: fill variables with zeros.")


    (options, args) = parser.parse_args()
    
    if len(args) == 4 :
        stnFile, g5File, iso_t1, iso_t2 = args
        t1, t2 = (isoparser(iso_t1), isoparser(iso_t2))
    elif len(args) == 2 :
        stnFile, g5File = args
        t1, t2 = (None,None)
    else:
        parser.error("must have 2 or 4 arguments: stnFile g5File [iso_t1 iso_t2]")

    if options.Vars is not None:
        options.Vars = options.Vars.split(',')

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
            
    # Station Locations
    # -----------------
    stnName, stnLon, stnLat = getStations(stnFile)

    # Time range
    # ----------
    Dt = timedelta(seconds=options.dt_secs)
    tyme = getTyme(Dt,t1,t2)

    # Get Variables and Metadata
    # --------------------------
    Vars, levs, levUnits = getVars(g5File)
    
    # Write output file
    # -----------------
    writeNC(stnName,stnLon,stnLat,tyme,Vars,levs,levUnits,options)

