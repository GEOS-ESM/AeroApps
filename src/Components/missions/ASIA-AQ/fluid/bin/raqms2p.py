#!/usr/bin/env python
#
# Convert RAQMS output to COARDS compliant file on pressure coordinates.
# The input grid is already regular lat/lon.
#
# Arlindo da Silva, March 2019.

from dateutil.parser import parse as isoparser
from optparse   import OptionParser
from datetime   import datetime, timedelta
from netCDF4    import Dataset
import numpy as npy
import os

from curv2llp   import *

#---
def getTyme(ncin):
    t = ncin.variables['time']
    if len(t) != 1:
        raise ValueError, 'Fix me, I can only handle 1 step per file for now'
    yyyymmddhh = t.units
    tyme = datetime(int(yyyymmddhh[0:4]),
                    int(yyyymmddhh[4:6]),
                    int(yyyymmddhh[6:8]),
                    int(yyyymmddhh[8:10]),
                    )
    return npy.array([tyme,])

#---
def get_cP(n,ncin):
    """
    Given time level "n", return 3D pressure for this time.
    """
    delp = nc.variables['delp'][n,:,:,:]
    nz, ny, nx = delp.shape
    
    pe = npy.zeros((nz+1,ny,nx))
    pe[0,:,:] = 0.4 # derived from comparison with ps.
    for k in range(nz):
        pe[k+1,:,:] = pe[k,:,:] + delp[k,:,:]

    cP = delp # to conserve memory
    for k in range(nz):
        cP[k,:,:] = (pe[k,:,:]+pe[k+1,:,:])/2.

    return cP


#---
class myCurv2LLP(Curv2LLP):
    
    def get_cP(self,n,ncin):
        """
        Given time level "n", return 3D pressure for this time.
        """
        delp = nc.variables['delp'][n,:,:,:]
        nz, ny, nx = delp.shape
        
        pe = npy.zeros((nz+1,ny,nx))
        pe[0,:,:] = 0.4 # derived from comparison with ps.
        for k in range(nz):
            pe[k+1,:,:] = pe[k,:,:] + delp[k,:,:]

        cP = delp # to conserve memory
        for k in range(nz):
            cP[k,:,:] = (pe[k,:,:]+pe[k+1,:,:])/2.

        return cP

def append_vars(ncin, zlib=False):
    """
    Append AOD, total column PM2.5, OC, CO to original file.
    """

    # Need this for getting AOD
    # -------------------------
    delp = 100*ncin.variables['delp'][0,:,:,:] # in Pa
    g = 9.80665

    # AOD 
    # ---
    rho = ncin.variables['rho'][0,:,:,:]
    for v in ['ext_tot', 'ext_sulf', 'ext_bcoc', 'ext_dust', 'ext_salt']:
        s = v.replace('ext_','')
        aod = npy.sum(ncin.variables[v][0,:,:,:]*delp[:,:,:]/rho[:,:,:],axis=0)/g
        this = ncin.createVariable('aod_'+s,'f4',('time', 'south_north', 'west_east'),zlib=zlib)
        this.long_name = u'%s Aerosol Optical Depth'%s.capitalize()
        this.missing_value = UNDEF # out of domain
        this.units = '1'
        this[0,:,:] = aod[:,:]

    # PM2.5
    # -----
    oc1 = ncin.variables['ioc1'][0,:,:,:] 
    so4aer = ncin.variables['iso4aer'][0,:,:,:] 
    oc1 = ncin.variables['ioc1'][0,:,:,:] 
    oc2 = ncin.variables['ioc2'][0,:,:,:] 
    bc1 = ncin.variables['ibc1'][0,:,:,:] 
    bc2 = ncin.variables['ibc2'][0,:,:,:] 
    du1 = ncin.variables['idu1'][0,:,:,:] 
    du2 = ncin.variables['idu2'][0,:,:,:] 
    ss1 = ncin.variables['iss1'][0,:,:,:] 
    ss2 = ncin.variables['iss2'][0,:,:,:] 

    bcoc=bc1+bc2+oc1+oc2
    ssalt=ss1+0.942*ss2
    dust=du1+0.286*du2

    # Assume all sulfate is ammonium sulfate
    amo = 132.
    amd = 28.9644
    vmmr = amo/amd
    so4aer=1e9*rho*so4aer*vmmr
    # carbon aerosols
    amo = 16.8
    amd = 28.9644
    vmmr = amo/amd
    bcoc=1e9*rho*bcoc*vmmr
    # dust
    dust=1e9*rho*dust
    # ssalt
    ssalt=1e9*rho*ssalt

    pm25=so4aer+bcoc+ssalt+dust
    this = ncin.createVariable('pm25','f4',('time','top_bottom', 'south_north', 'west_east'),zlib=zlib)
    this.long_name = 'pm2.5'
    this.missing_value = UNDEF # out of domain
    this.units = 'ug m**-3'
    this[0,:,:,:] = pm25[:,:,:]

 
    # total cloumn PM2.5
    # ------------------
    nz, ny, nx = pm25.shape
    var = npy.zeros((ny,nx))
    var = npy.sum(pm25[:,:,:]/1e9*delp[:,:,:]/g/rho, axis=0)
    this = ncin.createVariable('tcpm25','f4',('time', 'south_north', 'west_east'),zlib=zlib)
    this.long_name = 'total column pm2.5'
    this.missing_value = UNDEF # out of domain
    this.units = 'kg m**-2'
    this[0,:,:] = var[:,:]

    # total cloumn OC and CO
    # ----------------------
    co = ncin.variables['ico'][0,:,:,:]*28.0/28.97  # kg/kg
    var = npy.zeros((ny,nx))
    var = npy.sum(co[:,:,:]*delp[:,:,:]/g, axis=0)
    this = ncin.createVariable('tcco','f4',('time', 'south_north', 'west_east'),zlib=zlib)
    this.long_name = 'total column co'
    this.missing_value = UNDEF # out of domain
    this.units = 'kg m**-2'
    this[0,:,:] = var[:,:]
    
    oc =( ncin.variables['ioc1'][0,:,:,:] + ncin.variables['ioc2'][0,:,:,:] )*12.0/28.97  # kg/kg
    var = npy.zeros((ny,nx))
    var = npy.sum(oc[:,:,:]*delp[:,:,:]/g, axis=0)
    this = ncin.createVariable('tcoc','f4',('time', 'south_north', 'west_east'),zlib=zlib)
    this.long_name = 'total column oc'
    this.missing_value = UNDEF # out of domain
    this.units = 'kg m**-2'
    this[0,:,:] = var[:,:]



def append_lml(ncin, options, zlib=False):
    """
    Append vars at lowest model level to netcdf file.
    """

    # Open output NC file
    # -------------------
    nc = Dataset(options.outFile,'a',format=options.format)

    for v in ['iso4aer','ioc1','ioc2','ibc1','ibc2','idu1','idu2','idu3','idu4','idu5',
        'iss1','iss2','iss3','iss4','iss5','ico','no','ino2','iso2','o3vmr','rho','pm25']:
        this = nc.createVariable(v+'_lml','f4',('time','lat','lon'),zlib=zlib)
        this.long_name = ''
        this.missing_value = UNDEF # out of domain
        this.units = ''
        this[0,:,:] = ncin.variables[v][0,34,:,:]

    nc.close()    

def append_density(ncin, zlib=False):
  
    # append air density to input file 

    cP = get_cP(0,ncin)
    psfc = ncin.variables['psfc'][0,:,:] # in hPa 
    a = 0.286
    theta = ncin.variables['ttheta'][0,:,:,:]
    nz, ny, nx = theta.shape
    rho = theta # to conserve memory
    Rd = 287.05 
    for k in range(nz):
        T = theta[k,:,:]/pow((psfc/cP[k,:,:]),a) 
        rho[k,:,:] = cP[k,:,:]*100.0/Rd/T

    this = ncin.createVariable('rho','f4',('time', 'top_bottom', 'south_north', 'west_east'),zlib=zlib)
    this.long_name = 'air density'
    this.missing_value = UNDEF # out of domain
    this.units = 'kg m-3'
    this[0,:,:,:] = rho[:,:,:]
    
        

#---        
if __name__ == "__main__":

    format = 'NETCDF4'
    algo = 'linear'
#   rP = '925,850,700,600'
    rP = '1000,975,950,925,900,875,850,825,800,775,750,725,700,650,600,550,500,450,400,350,300,250,200,150,100'

#   Parse command line options
#   --------------------------
    parser = OptionParser(usage="Usage: %prog [OPTIONS] raqms_File",
                          version='1.0.0' )

    parser.add_option("-o", "--output", dest="outFile", default=None,
              help="Output NetCDF file (default: same was input with uwhyb replaced with uwllp)")

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
        raqms_File = args[0]
    else:
        parser.error("must have 1 argument: raqms_File")

    if options.Vars is not None:
        options.Vars = options.Vars.split(',')

    if options.outFile is None:
        options.outFile = raqms_File.replace('uwhyb_','uwllp_')
    if options.outFile == raqms_File:
        options.outFile = 'raqms_llp.nc4' # do not overwite input
        
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
    
    interm = name + '.tmp.nc'
    os.system("cp "+raqms_File+" "+interm)
    print "cp "+raqms_File+" "+interm


    # Open the input file
    # -------------------
   # nc = Dataset(raqms_File)
    nc = Dataset(interm,'a',format=options.format)

    # Time range
    # ----------
    tyme = getTyme(nc)
    
    # calculate air density 
    # ---------------------
    append_density(nc, zlib=False)

    # append pm2.5, aod, total column OC and CO 
    # -----------------------------------------
    append_vars(nc, zlib=False)

    # Coordinates
    # -----------
    lon = nc.variables['lon'][:]
    lat = nc.variables['lat'][:]

    # Instantiate regridding class
    # ----------------------------
    r = myCurv2LLP(lon,lat,options.rP)
    
    # Write output file
    # -----------------
    r.writeNC ( tyme, options, nc, zlib=False )

    # append variables at the lowest model level 
    append_lml(nc, options, zlib=False)

    os.remove(interm)
