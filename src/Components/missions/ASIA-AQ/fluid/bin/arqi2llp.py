#!/usr/bin/env python
#
# Convert ARQI output to COARDS compliant file on a regular lat-lon-pressure grid.
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
        lon = nc.variables['rlon1']
        lat = nc.variables['rlat1']
        pref = nc.variables['pref']
        a = ncin.variables['a_1'][:]
        b = ncin.variables['b_1'][:]
        nx, ny, nz = lon.shape[0], lat.shape[0], len(a)
        ps = nc.variables['P0'][n,:,:]*100.0  # Pa
        cP = npy.zeros((nz,ny,nx))
        for k in range(nz):
            cP[k,:,:] = npy.exp(a[k] + b[k]*npy.log(ps[:,:]/pref))/100.  #hPa 
        return cP
        

#---
def getTyme(ncin):
    T = ncin.variables['time']
    # hourss since 2018-08-07 00:00:00
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
def append_vars(ncin, options, zlib=False):
    """
    Append total column PM2.5, OC, CO 
    """
    g = 9.80665
    nc = Dataset(options.outFile,'a',format=options.format)
    p = myCurv2LLP.get_cP(0,ncin)
    nz, ny, nx = p.shape
    dp= npy.zeros((nz,ny,nx))


#---
def append_lml( nc, zlib=False):
    """
    Append vars at lowest model level to the oroginal file.
    """

    for v in ['AF','AC','AMAC', 'AM25','NIAC', 'NI25','SUAC', 'SU25', 
        'PCAC', 'PC25','OCAC', 'OC25','ECAC', 'EC25','CMAC', 'CM25',
        'SSC5', 'SSC6', 'SSC7', 'SSC8', 'SSC9', 'SSCA', 'SSCB',
        'C3H8','CO','HCHO','ISOP','N2','NH3','NO','O3','OH','S2','PAN']:

        print 'Appending ', v
        this = nc.createVariable(v+'_lml','f4',('time','rlat1', 'rlon1'),zlib=zlib)
        this.long_name = ''
        this.missing_value = UNDEF # out of domain
        this.units = ''

        nz, ny, nx = nc.variables[v][0,:,:,:].shape
        this[0,:,:] = nc.variables[v][0,nz-1,:,:] #  top to bottom 


   
#---        
if __name__ == "__main__":

    format = 'NETCDF4'
    algo = 'linear'
    rP = '925,850,700,600'

#   Parse command line options
#   --------------------------
    parser = OptionParser(usage="Usage: %prog [OPTIONS] arqi_File",
                          version='1.0.0' )

    parser.add_option("-o", "--output", dest="outFile", default=None,
              help="Output NetCDF file (default: same was input with arqiout replaced with arqillp)")

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
        arqi_File = args[0]
    else:
        parser.error("must have 1 argument: arqi_File")

    if options.Vars is not None:
        options.Vars = options.Vars.split(',')

    if options.outFile is None:
        dirn = os.path.dirname(arqi_File)
        base = os.path.basename(arqi_File)
        options.outFile = dirn + '/arqi.llp.' + base.replace('_uconversion_particles.netcdf4.compressed','.nc4')
    if options.outFile == arqi_File:
        options.outFile = 'arqi_llp.nc4' # do not overwite input
        
    options.rP = [ float(p) for p in options.rP.split(',') ]
    
    if options.ignore is not None:
        options.ignore = options.ignore.split(',')
    else:
        options.ignore = []
        
    options.ignore += ['rlon1', 'rlat1', 'pref', 'reftime', 'rotate_pole',
                       'a', 'b','lon_1', 'lat_1', 'p', 'time'] # these are coordinates

    # Create consistent file name extension
    # -------------------------------------
    name, ext = os.path.splitext(options.outFile)
    if 'NETCDF4' in options.format:
        options.outFile = name + '.nc4'
    elif 'NETCDF3' in options.format:
        options.outFile = name + '.nc'
    else:
        raise ValueError, 'invalid extension <%s>'%ext

    # copy original file 
    # ------------------

    tpath = os.path.dirname(name)
    tname = os.path.basename(arqi_File)
    interm = tpath + os.sep + tname + '.tmp.nc'
    
    os.system("cp "+arqi_File+" " + interm)

    # Open the input file
    # -------------------
    nc = Dataset(interm,'a',format=options.format)

    # Time range
    # ----------
    tyme = getTyme(nc)

    # Append lml vars
    # ---------------
    append_lml(nc, zlib=False)

    
    # Instantiate regridding class
    # ----------------------------
    cLon = nc.variables['lon_1'][:,:]
    cLat = nc.variables['lat_1'][:,:]
    r = myCurv2LLP(cLon,cLat,options.rP)
    
    # Write output file
    # -----------------
    r.writeNC ( tyme, options, nc, zlib=False )

    os.remove(interm)
