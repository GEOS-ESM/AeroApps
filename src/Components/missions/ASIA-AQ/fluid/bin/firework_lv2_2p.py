#!/usr/bin/env python
#
# Convert FireWork output to COARDS compliant file on a regular lat-lon-pressure grid.
#
# Arlindo da Silva, March 2019.

from dateutil.parser import parse as isoparser
from optparse   import OptionParser
from datetime   import datetime, timedelta
from netCDF4    import Dataset
import numpy as npy
import os

from curv2llp   import *

class myCurv2LLP(Curv2LLP):
    
    def get_cP(self,n,ncin):
        """
        Given time level "n", return 3D pressure for this time.
        """
        lon = nc.variables['lon']
        lat = nc.variables['lat']
        pref = nc.variables['pref']
        a = ncin.variables['a_2'][:]
        b = ncin.variables['b_2'][:]
        nx, ny, nz = lon.shape[0], lat.shape[0], len(a)
        ps = nc.variables['P0'][n,:,:]*100.0  # Pa
        cP = npy.zeros((nz,ny,nx))
        for k in range(nz):
            cP[k,:,:] = npy.exp(a[k] + b[k]*npy.log(ps[:,:]/pref))/100.  #hPa 
        print cP.min() 
        return cP
        

#---
def getTyme(ncin):
    T = ncin.variables['time']
    # e.g. hourss since 2018-08-07 00:00:00
    units = T.units
    isoT0 = units.split(' since ')[1].replace(' ','T')
    cDT = units.split(' ')[0]
    sec = timedelta(seconds=1)
    if cDT == "hours":
        DT = 60
    else:
        raise ValueError, "fix me, can only handle hours as unit of time"
    t0 = isoparser(isoT0)
    tyme = npy.array([t0 + int(t * DT)*sec for t in T])
    return tyme

   
#---        
if __name__ == "__main__":

    format = 'NETCDF4'
    algo = 'linear'
    # rP = '925,850,700,600'
    rP = '1000,975,950,925,900,875,850,825,800,775,750,725,700,650,600,550,500,400,300,200'

#   Parse command line options
#   --------------------------
    parser = OptionParser(usage="Usage: %prog [OPTIONS] firew_File",
                          version='1.0.0' )

    parser.add_option("-o", "--output", dest="outFile", default=None,
              help="Output NetCDF file (default: same was input with firewout replaced with firewllp)")

    parser.add_option("-a", "--algorithm", dest="algo", default=algo,
              help="Interpolation algorithm, one of linear, cubic (default=%s)"\
                          %algo)

    parser.add_option("-V", "--vars", dest="Vars", default=None,
              help="Variables to sample (default=All)")

    parser.add_option("-x", "--exclude", dest="ignore", default=None,
              help="Variables to ignore (default=None)")
    
    parser.add_option("-f", "--format", dest="format", default=format,
              help="Output file format: one of NETCDF4, NETCDF4_CLASSIC, NETCDF3_CLASSIC or NETCDF3_64BIT (default=%s)"%format )

    parser.add_option("-p", "--levels", dest="rP", default=rP,
              help="Levels to sample (default=%s)"%rP)

    parser.add_option("-v", "--verbose",
                      action="store_true", dest="verbose",
                      help="Verbose mode.")

    parser.add_option("-n", "--dryrun",
                      action="store_true", dest="dryrun",
                      help="Dry-run mode: fill variables with zeros.")
    

    (options, args) = parser.parse_args()
    
    if len(args) == 1 :
        firew_File = args[0]
    else:
        parser.error("must have 1 argument: firew_File")

    if options.Vars is not None:
        options.Vars = options.Vars.split(',')

    if options.outFile is None:
        dirn = os.path.dirname(firew_File)
        base = os.path.basename(firew_File)
        options.outFile = dirn + '/firew.llp.' + base.replace('_uconversion_particles.netcdf4.compressed','.nc4')
    if options.outFile == firew_File:
        options.outFile = 'firew_llp.nc4' # do not overwite input
        
    options.rP = [ float(p) for p in options.rP.split(',') ]
    
    if options.ignore is not None:
        options.ignore = options.ignore.split(',')
    else:
        options.ignore = []
        
    options.ignore += ['lon', 'lat', 'pref',
                       'a_1', 'b_1','lon_1', 'lat_1', 'level1',
                       'a_2', 'b_2','level2', 
                       'a_3', 'b_3','level3',
                       'a_5', 'b_5','level5',
                       'time', 'GZ_vgrid5'] 

    options.Vars = ['NH3','C3H8','HCHO','ISOP','CO', 'OH', 'PAN','RHO']

    # Create consistent file name extension
    # -------------------------------------
    name, ext = os.path.splitext(options.outFile)
    if 'NETCDF4' in options.format:
        options.outFile = name + '.nc4'
    elif 'NETCDF3' in options.format:
        options.outFile = name + '.nc'
    else:
        raise ValueError, 'invalid extension <%s>'%ext


    # Open the input file
    # -------------------
    nc = Dataset(firew_File,format=options.format)

    # Time range
    # ----------
    tyme = getTyme(nc)

    
    # Instantiate regridding class
    # ----------------------------
    cLon = nc.variables['lon'][:]
    cLat = nc.variables['lat'][:]
    r = myCurv2LLP(cLon,cLat,options.rP)
    
    # Write output file
    # -----------------
    r.writeNC ( tyme, options, nc, zlib=False )

    
