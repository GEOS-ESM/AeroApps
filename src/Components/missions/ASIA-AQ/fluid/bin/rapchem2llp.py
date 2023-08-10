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
        press = (ncin.variables['p'][n,:,:,:])/100. # in hPa
        return press[:,:,:]

#---
def getTyme(ncin):
    t = ncin.variables['Times']

    if len(t[:,0]) != 1:
        raise ValueError, 'Fix me, I can only handle 1 step per file for now'
    tstr=''.join(t[0,:])
    print(tstr[0:4])   # 2018-08-07_06:00:00
    tyme = datetime(int(tstr[0:4]),
                    int(tstr[5:7]),
                    int(tstr[8:10]),
                    int(tstr[11:13])) 
    return npy.array([tyme,])
    
#---
def append_lml(nc, zlib=False):
    """
    Append vars at lowest model level to the oroginal file.
    """
    for v in ['PM2_5_DRY', 'nh4aj', 'nh4ai','no3aj', 'no3ai','sulf', 'so4aj', 'so4ai',
        'orgpai', 'orgpaj','ecj', 'eci','seas','asoa1j', 'asoa1i','asoa2j', 'asoa2i','asoa3j', 
        'asoa3i','asoa4j', 'asoa4i','bsoa1j', 'bsoa1i','bsoa2j', 'bsoa2i','bsoa3j', 'bsoa3i', 
        'bsoa4j', 'bsoa4i','soila','smoke','co','no2','no','o3','so2','ALT']:
        var = nc.variables[v]
        this = nc.createVariable(v+'_lml','f4',('Time', 'south_north', 'west_east'),zlib=zlib)
        this.description = var.description+' at lowest model level'
        this.missing_value = UNDEF # out of domain
        this.units = var.units
        this[0,:,:] = var[0,0,:,:] 

#---
def append_vars(nc, zlib=False): 
    """
    Append total column vars 
    """
   
    co = nc.variables['co'][0,:,:,:]  # ppmv 
    co = co/1e6*28.0/28.97 # kg kg-1
    nz, ny, nx = co.shape
    oc = npy.zeros((nz,ny,nx))
    for v in ['asoa1j', 'asoa1i','asoa2j', 'asoa2i','asoa3j',
        'asoa3i','asoa4j', 'asoa4i','bsoa1j', 'bsoa1i','bsoa2j', 'bsoa2i','bsoa3j', 'bsoa3i',
        'bsoa4j', 'bsoa4i','orgpai', 'orgpaj']:
        oc = oc + nc.variables[v][0,:,:,:]  # ug kg-1 
    rho = 1.0/nc.variables['ALT'][0,:,:,:]  # kg m-3 
    dz = nc.variables['dz'][0,:,:,:]  # m 
    pm25 = nc.variables['PM2_5_DRY'][0,:,:,:]  # ug m-3 

    # total cloumn OC
    #----------------
    this = nc.createVariable('tcoc','f4',('Time', 'south_north', 'west_east'),zlib=zlib)
    this.description = 'total column OC'
    this.missing_value = UNDEF # out of domain
    this.units = 'kg m-2'
    this[0,:,:] = npy.sum(oc*rho*dz/1e9, axis=0) # kg m-2 

    # total cloumn CO 
    #----------------
    this = nc.createVariable('tcco','f4',('Time', 'south_north', 'west_east'),zlib=zlib)
    this.description = 'total column CO'
    this.missing_value = UNDEF # out of domain
    this.units = 'kg m-2'
    this[0,:,:] = npy.sum(co*rho*dz, axis=0)  # kg m-2

    # total cloumn PM2.5
    #----------------
    this = nc.createVariable('tcpm25','f4',('Time', 'south_north', 'west_east'),zlib=zlib)
    this.description = 'total column CO'
    this.missing_value = UNDEF # out of domain
    this.units = 'kg m-2'
    this[0,:,:] = npy.sum(pm25*dz/1e9, axis=0)  # kg m-2


#---        
if __name__ == "__main__":

    format = 'NETCDF4'
    algo = 'linear'
    rP = '925,850,700,600'

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
        
    options.ignore += ['XLONG', 'XLAT', 'p', 'time','lon','lat'] # these are coordinates

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
    os.system("cp "+wrf_File+" interm.nc")

    # Open the input file
    # -------------------
    nc = Dataset( "interm.nc", 'a', format=options.format)

    # Append vars 
    # -----------
    append_lml(nc,zlib=False)
    append_vars(nc,zlib=False)

    # Time range
    # ----------
    tyme = getTyme(nc)
    
    
    # Instantiate regridding class
    # ----------------------------
    cLon = npy.squeeze(nc.variables['XLONG'][:,:,:])
    cLat = npy.squeeze(nc.variables['XLAT'][:,:,:])
    r = myCurv2LLP(cLon,cLat,options.rP)
    
    # Write output file
    # -----------------
    r.writeNC ( tyme, options, nc, zlib=False )
    
