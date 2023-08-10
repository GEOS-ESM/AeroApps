#!/usr/bin/env python
#
# Append 3D PM2.5, total column PM2.5, total column OC 
# can only handle 1 time step per file for now


from dateutil.parser import parse as isoparser
from optparse   import OptionParser
from datetime   import datetime, timedelta
from netCDF4    import Dataset
import numpy as npy
import os

def append_vars(ncin, options, zlib=False):
    """
    Append 3-d PM2.5, total column PM2.5, OC
    """

    g = 9.80665
    
    nc = Dataset(options.outFile,'a',format=options.format)
    rho = ncin.variables['p3089'][0,:,:,:]  # air density in kg m-3 
    ps = npy.exp(ncin.variables['lnsp'][0,:,:])  # surface pressure  
    p = ncin.variables['level'][:] # pressure levels, top to bottom in Pa 
    nz, ny, nx = rho.shape
    dp = npy.zeros((nz,ny,nx))
    dp[0,:,:]=(p[0]+p[1])/2. 
    for k in range(1,nz-1):
        dp[k,:,:] = (p[k+1]-p[k-1])/2. 
    dp[nz-1,:,:] = ps - (p[nz-1]+p[nz-2])/2.
   

    # PM2.5 = RHO * ( 1 * SS1 / 4.3 + 0.5 * SS2 / 4.3 + 1 * DD1 + 1 * DD2 + 0.7 * OM1 + 0.7 * OM2 + O.7 * SU1 + 1 * BC1 + 1 * BC2 )
    ss1 = ncin.variables['aermr01'][0,:,:,:]  
    ss2 = ncin.variables['aermr02'][0,:,:,:]  
    dd1 = ncin.variables['aermr04'][0,:,:,:]  
    dd2 = ncin.variables['aermr05'][0,:,:,:]  
    om1 = ncin.variables['aermr07'][0,:,:,:]  
    om2 = ncin.variables['aermr08'][0,:,:,:]  
    su1 = ncin.variables['aermr11'][0,:,:,:]  
    bc1 = ncin.variables['aermr09'][0,:,:,:]  
    bc2 = ncin.variables['aermr10'][0,:,:,:]  
    pm2p5_3d = rho*(ss1/4.3 + 0.5*ss2/4.3 + dd1 + dd2 + 0.7*om1 + 0.7*om2 + 0.7*su1 + bc1 + bc2 )
    tcpm2p5 = npy.sum( pm2p5_3d[:,:,:]*dp/rho/g, axis=0)
    tcom = npy.sum( (om1[:,:,:]+om2[:,:,:])*dp/g, axis=0)
    
    this = nc.createVariable('pm2p5_3d','f4',('time', 'level', 'latitude', 'longitude'),zlib=zlib)
    this.long_name = 'pm 2.5'
    # this.missing_value = UNDEF # out of domain
    this.units = 'kg m**-3'
    this[0,:,:,:] = pm2p5_3d[:,:,:]
    
    this = nc.createVariable('tcpm2p5','f4',('time','latitude', 'longitude'),zlib=zlib)
    this.long_name = 'total column pm 2.5'
    # this.missing_value = UNDEF # out of domain
    this.units = 'kg m**-2'
    this[0,:,:] = tcpm2p5[:,:]

    this = nc.createVariable('tcom','f4',('time','latitude', 'longitude'),zlib=zlib)
    this.long_name = 'total column organic matter'
    # this.missing_value = UNDEF # out of domain
    this.units = 'kg m**-2'
    this[0,:,:] = tcom[:,:]

    nc.close()    

#---        
if __name__ == "__main__":

    format = 'NETCDF4'

#   Parse command line options
#   --------------------------
    parser = OptionParser(usage="Usage: %prog [OPTIONS] cams_File",
                          version='1.0.0' )

    parser.add_option("-o", "--output", dest="outFile", default=None,
              help="Output NetCDF file (default: same as input but replace merged with prep)")

    parser.add_option("-V", "--vars", dest="Vars", default=None,
              help="Variables to sample (default=All)")

    parser.add_option("-x", "--exclude", dest="ignore", default=None,
              help="Variables to ignore (default=None)")
    
    parser.add_option("-f", "--format", dest="format", default=format,
              help="Output file format: one of NETCDF4, NETCDF4_CLASSIC, NETCDF3_CLASSIC or NETCDF3_64BIT (default=%s)"%format )

    parser.add_option("-v", "--verbose",
                      action="store_true", dest="verbose",
                      help="Verbose mode.")

    parser.add_option("-n", "--dryrun",
                      action="store_true", dest="dryrun",
                      help="Dry-run mode: fill variables with zeros.")


    (options, args) = parser.parse_args()
    
    if len(args) == 1 :
        cams_File = args[0]
    else:
        parser.error("must have 1 argument: cams_File")

    if options.Vars is not None:
        options.Vars = options.Vars.split(',')

    if options.outFile is None:
        options.outFile = cams_File.replace('merged_','prep_')
        
    

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
    os.system("cp "+cams_File+" "+options.outFile)

    # Open the input file
    # -------------------
    nc = Dataset(cams_File)

    append_vars(nc, options, zlib=False)



    
