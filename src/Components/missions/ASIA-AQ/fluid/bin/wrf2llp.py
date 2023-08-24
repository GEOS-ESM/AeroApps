#!/usr/bin/env python
#
# Convert WRFchem output to COARDS compliant file on a regular lat-lon-pressure grid.
#
# Arlindo da Silva, March 2019.

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
        return ncin.variables['p'][n,:,:,:]/100.  # in hPa

#---
def getTyme(ncin):
    t = ncin.variables['time']
    units = t.units
    hh = units.split('Z')[0]
    yyyymmdd = units.split(' ')[1]
    tyme = datetime(int(yyyymmdd[0:4]),
                    int(yyyymmdd[4:6]),
                    int(yyyymmdd[6:8]),
                    int(hh))
    return npy.array([tyme,])
    
#---
def append_lml(nc, zlib=False):
    """
    Append vars at lowest model level to the oroginal file.
    """
    for v in ['oc', 'bc', 'pm25']:
        var = nc.variables[v]
        this = nc.createVariable(v+'_lml','f4',('time', 'south_north', 'west_east'),zlib=zlib)
        this.long_name = var.long_name+' at lowest model level'
        this.missing_value = UNDEF # out of domain
        this.units = var.units
        this[0,:,:] = var[0,0,:,:] 
        
#---
def append_vars(nc, zlib=False):
    """
    Append total column vars 
    """
    # total cloumn OC
    #----------------
    dz = nc.variables['dz'][0,:,:,:]  # in m 
    nz, ny, nx = dz.shape
    tcoc = npy.zeros((ny, nx))
    tcoc = npy.sum(nc.variables['oc'][0,:,:,:]*dz/1e9, axis=0)  # in kg m-2
    this = nc.createVariable('tcoc','f4',('time', 'south_north', 'west_east'),zlib=zlib)
    this.long_name = 'total column OC'
    this.missing_value = UNDEF # out of domain
    this.units = 'kg m-2'
    this[0,:,:] = tcoc[:,:] 
    

    # total cloumn pm2.5
    #----------------
    tc = npy.zeros((ny,nx))
    tc = npy.sum(nc.variables['pm25'][0,:,:,:]*dz/1e9, axis=0 ) 
    this = nc.createVariable('tcpm25','f4',('time', 'south_north', 'west_east'),zlib=zlib)
    this.long_name = 'total column PM2.5'
    this.missing_value = UNDEF # out of domain
    this.units = 'kg m-2'
    this[0,:,:] = tc[:,:] 


#---        
if __name__ == "__main__":

    format = 'NETCDF4'
    algo = 'linear'
    #rP = '925,850,700,600'
    rP = '1000,975,950,925,900,875,850,825,800,775,750,725,700,650,600,550,500'
    #rP = '1000,975,950,925,900,875,850,825,800,775,750,725,700,650,600,550,500,450,400,350,300,250,200,150,100'

#   Parse command line options
#   --------------------------
    parser = OptionParser(usage="Usage: %prog [OPTIONS] wrf_File",
                          version='1.0.0' )

    parser.add_option("-o", "--output", dest="outFile", default=None,
              help="Output NetCDF file (default: same was input with wrfout replaced with wrfllp)")

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
        wrf_File = args[0]
    else:
        parser.error("must have 1 argument: wrf_File")

    if options.Vars is not None:
        options.Vars = options.Vars.split(',')

    if options.outFile is None:
        options.outFile = wrf_File.replace('wrfout_','wrfllp_')
    if options.outFile == wrf_File:
        options.outFile = 'wrfllp.nc4' # do not overwite input
        
    options.rP = [ float(p) for p in options.rP.split(',') ]
    
    if options.ignore is not None:
        options.ignore = options.ignore.split(',')
    else:
        options.ignore = []
        
    options.ignore += ['lon', 'lat', 'p', 'time'] # these are coordinates

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
    os.system("cp "+wrf_File+" "+wrf_File+".interm.nc")

    # Open the input file
    # -------------------
    nc = Dataset(wrf_File+".interm.nc",'a',format=options.format)

    # Time range
    # ----------
    tyme = getTyme(nc)

    # Append vars 
    # -----------
    print("append vars ... ")
    append_lml(nc,zlib=False)
    append_vars(nc,zlib=False)
    print("finished append vars ... ")
    

    # Instantiate regridding class
    # ----------------------------
    print("regridding ... ")
    cLon = npy.squeeze(nc.variables['lon'][:,:,:])
    cLat = npy.squeeze(nc.variables['lat'][:,:,:])
    r = myCurv2LLP(cLon,cLat,options.rP)
    
    # Write output file
    # -----------------
    print("write output ... ")
    r.writeNC ( tyme, options, nc, zlib=False )
    os.system("rm "+wrf_File+".interm.nc")


