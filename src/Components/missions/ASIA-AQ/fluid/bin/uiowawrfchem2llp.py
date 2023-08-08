#!/usr/bin/env python
#
# Convert UIOWA WRFchem output to COARDS compliant file on a regular lat-lon-pressure grid.
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
        return ncin.variables['P'][0,:,:,:]/100.  # in hPa

#---
def getTyme(wrf_File):
    # input  e.g.:  WRFCHEM_UIOWA_2019-07-22_12.nc 
    base = os.path.basename(wrf_File)
    tyme = datetime(int(base[14:18]),
                    int(base[19:21]),
                    int(base[22:24]),
                    int(base[25:27]))
    return npy.array([tyme,])
    
#---
def append_lml(nc, zlib=False):
    """
    Append vars at lowest model level to the oroginal file.
    """
    for v in ['TEMP', 'RHO', 'RH', 'CLDFRA', 'U', 'V', 
        'QVAPOR', 'CO', 'HCHO', 'NO', 'NO2', 'NH3', 'O3', 'SO2', 
        'PM25', 'NITRATE', 'OC', 'BC', 'SULF', 'OIN', 'AMMONIUM', 'NA', 'CL']:
        var = nc.variables[v]
        this = nc.createVariable(v+'_lml','f4',('Time', 'lat', 'lon'),zlib=zlib)
        this.description = v +' at lowest model level'
        this.missing_value = UNDEF # out of domain
        if 'units' in var.ncattrs():
            this.units = var.units
        this[0,:,:] = var[0,0,:,:] 
        
#---
def append_vars(nc, zlib=False):
    """
    Append total column vars 
    """
    g = 9.80665
    rho = nc.variables['RHO'][0,:,:,:]  # kg m-3 
    psfc = nc.variables['PSFC'][0,:,:]  # Pa 
    p = nc.variables['P'][0,:,:,:]  # Pa 
    nz, ny, nx = p.shape 
    dz = npy.zeros((nz,ny,nx))
    for k in range(nz):   # bottom to top 
        if k == 0:
            dz[k,:,:]= (psfc-p[0,:,:])/rho[0,:,:]/g 
        else: 
            dz[k,:,:]= (p[k-1,:,:]-p[k,:,:])/rho[k,:,:]/g 

    # total cloumn OC
    #----------------
    tcoc = npy.zeros((ny, nx))
    tcoc = npy.sum(nc.variables['OC'][0,:,:,:]*dz/1e9, axis=0)  # in kg m-2
    this = nc.createVariable('tcoc','f4',('Time', 'lat', 'lon'),zlib=zlib)
    this.description = 'total column OC'
    this.missing_value = UNDEF # out of domain
    this.units = 'kg m-2'
    this[0,:,:] = tcoc[:,:] 
    

    # total cloumn pm2.5
    #----------------
    tc = npy.zeros((ny,nx))
    tc = npy.sum(nc.variables['PM25'][0,:,:,:]*dz/1e9, axis=0 ) 
    this = nc.createVariable('tcpm25','f4',('Time', 'lat', 'lon'),zlib=zlib)
    this.description = 'total column PM2.5'
    this.missing_value = UNDEF # out of domain
    this.units = 'kg m-2'
    this[0,:,:] = tc[:,:] 


    # total cloumn CO
    #----------------
    tc = npy.zeros((ny,nx))
    tc = npy.sum(nc.variables['CO'][0,:,:,:]*28.0/28.97*rho*dz/1e6, axis=0 ) 
    this = nc.createVariable('tcco','f4',('Time', 'lat', 'lon'),zlib=zlib)
    this.description = 'total column CO'
    this.missing_value = UNDEF # out of domain
    this.units = 'kg m-2'
    this[0,:,:] = tc[:,:] 

#---        
if __name__ == "__main__":

    format = 'NETCDF4'
    algo = 'linear'
    #rP = '925,850,700,600'
    rP = '1000,975,950,925,900,875,850,825,800,775,750,725,700,650,600,550,500,400,300,200,100'
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
        
    options.ignore += ['lon', 'lat', 'P', 'lev'] # these are coordinates

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
    tyme = getTyme(wrf_File)
  

    # Append vars 
    # -----------
    print("append vars ... ")
    append_lml(nc,zlib=False)
    append_vars(nc,zlib=False)
    

    # Instantiate regridding class
    # ----------------------------
    print("regridding ... ")
    cLon = npy.squeeze(nc.variables['lon'][:,:])
    cLat = npy.squeeze(nc.variables['lat'][:,:])
    r = myCurv2LLP(cLon,cLat,options.rP)
    
    # Write output file
    # -----------------
    print("write output ... ")
    r.writeNC ( tyme, options, nc, zlib=False )
    os.system("rm "+wrf_File+".interm.nc")


