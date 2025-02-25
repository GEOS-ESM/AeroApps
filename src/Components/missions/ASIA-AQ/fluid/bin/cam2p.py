#!/usr/bin/env python
#
# Convert CAM output to COARDS compliant file on pressure coordinates.
# The input grid is already regular lat/lon.
#
# Arlindo da Silva, March 2019.

from dateutil.parser import parse as isoparser
from optparse   import OptionParser
from datetime   import datetime, timedelta
from netCDF4    import Dataset, num2date
import numpy as npy
import os

from curv2llp   import *

def get_cP(n,ncin):
    """
    Given time level "n", return 3D pressure for this time.
    """
    lon = nc.variables['lon']
    lat = nc.variables['lat']
    a = ncin.variables['hyam'][:]
    b = ncin.variables['hybm'][:]
    nx, ny, nz = lon.shape[0], lat.shape[0], len(a)
    ps = nc.variables['PS'][n,:,:]
    cP = npy.zeros((nz,ny,nx))
    for k in range(nz):
        cP[k,:,:] = (a[k] + b[k]*ps[:,:])/100.0  # hPa
    return cP

#---
class myCurv2LLP(Curv2LLP):
    
    def get_cP(self,n,ncin):
        """
        Given time level "n", return 3D pressure for this time.
        """
        lon = nc.variables['lon']
        lat = nc.variables['lat']
        a = ncin.variables['hyam'][:]
        b = ncin.variables['hybm'][:]
        nx, ny, nz = lon.shape[0], lat.shape[0], len(a)
        ps = nc.variables['PS'][n,:,:]
        cP = npy.zeros((nz,ny,nx))
        for k in range(nz):
            cP[k,:,:] = ( a[k] + b[k]*ps[:,:] )/100.0   # hPa
        return cP

def get_cPi(n,ncin):
    """
    Given time level "n", return 3D pressure (at model interfaces) for this time.
    """
    lon = nc.variables['lon']
    lat = nc.variables['lat']
    a = ncin.variables['hyai'][:]
    b = ncin.variables['hybi'][:]
    nx, ny, nz = lon.shape[0], lat.shape[0], len(a)
    ps = nc.variables['PS'][n,:,:]
    cPi = npy.zeros((nz,ny,nx))
    for k in range(nz):
        cPi[k,:,:] = (a[k] + b[k]*ps[:,:])/100.0  # hPa
    return cPi

#---
def getTyme(ncin):
    T = ncin.variables['time']
    # days since 2018-08-07 00:00:00
    # datetime = num2date(T[:],units=T.units, calendar=T.calendar)
    # print(datetime)
    units = T.units
    isoT0 = units.split(' since ')[1].replace(' ','T')
    cDT = units.split(' ')[0]
    sec = timedelta(seconds=1)
    if cDT == "days":
        DT = 24*60*60
    else:
        raise ValueError, "fix me, can only handle days as unit of time"
    t0 = isoparser(isoT0)
    tyme = npy.array([t0 + int(t * DT)*sec for t in T ])
    return tyme


#---
def append_density( ncin, zlib=False):
    this = ncin.createVariable('rho','f4',('time', 'lev', 'lat', 'lon'),zlib=zlib)
    this.long_name = 'air density'
    this.missing_value = UNDEF # out of domain
    this.units = 'kg m-3'
    tk = ncin.variables['T'][:,:,:,:]
    nt, nz, ny, nx = tk.shape
    
    for n in range(nt):
        cP = get_cP(n,ncin)
        rho = cP # to conserve memory
        Rd = 287.05 
        for k in range(nz):
            rho[k,:,:] = cP[k,:,:]*100.0/Rd/tk[n,k,:,:]
        this[n,:,:,:] = rho[:,:,:]


#---
def append_lml(ncin, options, zlib=False):
    """
    Append vars at lowest model level to netcdf file.
    """

    # Open output NC file
    # -------------------
    nc = Dataset(options.outFile,'a',format=options.format)

    for v in ['NH4','so4_a1','so4_a2','so4_a3','pom_a4','soa1_a1','soa1_a2','soa2_a1','soa2_a2',
       'soa3_a1','soa3_a2','soa4_a1','soa4_a2','soa5_a1','soa5_a2','bc_a1','bc_a4','dst_a1',
       'dst_a2','dst_a3','ncl_a1','ncl_a2','ncl_a3','CO','NO2','NO','O3','CO01','SO2','rho']:
        this = nc.createVariable(v+'_lml','f4',('time','lat','lon'),zlib=zlib)
        var = ncin.variables[v][:,:,:,:]
        this.long_name = ''
        this.missing_value = UNDEF # out of domain
        this.units = ''
        nt, nz, ny, nx = var.shapee
        for n in range(nt)
            this[n,:,:] = var[n,87,:,:]

    nc.close()

#---
def append_vars(ncin, options, zlib=False):
    """
    Append total column OC, CO to output file.
    """
    # Open output NC file
    # -------------------
    nc = Dataset(options.outFile,'a',format=options.format)

    g = 9.80665
    co = ncin.variables['CO'][:,:,:,:]*28.0/28.97 # kg kg-1
    nt, nz, ny, nx = co.shape
    oc = npy.zeros((nt,nz,ny,nx))

    for v in ['pom_a4','soa1_a1','soa1_a2','soa2_a1','soa2_a2',
        'soa3_a1','soa3_a2','soa4_a1','soa4_a2','soa5_a1','soa5_a2']: 
        oc = oc + ncin.variables[v][:,:,:,:]

    delp = npy.zeros((nt,nz,ny,nx))
    for n in range(nt):
        cPi = get_cPi(n,ncin)
        for k in range(nz):
            delp[n,k,:,:]= cPi[k+1,:,:]-cPi[k,:,:]  # hPa 


    # total column CO 
    # ---------------
    this = nc.createVariable('tcco','f4',('time', 'lat', 'lon'),zlib=zlib)
    this.long_name = 'total column CO'
    this.missing_value = UNDEF # out of domain
    this.units = 'kg m-2'
    for n in range(nt): 
        tcco = npy.sum(co[n,:,:,:]*delp[n,:,:,:]*100.0/g, axis=0)
        this[n,:,:] = tcco[:,:]

    # total column OC
    # ---------------
    this = nc.createVariable('tcoc','f4',('time', 'lat', 'lon'),zlib=zlib)
    this.long_name = 'total column OC'
    this.missing_value = UNDEF # out of domain
    this.units = 'kg m-2'
    for n in range(nt): 
        tcoc = npy.sum(oc[n,:,:,:]*delp[n,:,:,:]*100.0/g, axis=0)
        this[n,:,:] = tcoc[:,:]


#---        
if __name__ == "__main__":

    format = 'NETCDF4'
    algo = 'linear'
    rP = '925,850,700,600'

#   Parse command line options
#   --------------------------
    parser = OptionParser(usage="Usage: %prog [OPTIONS] cam_File",
                          version='1.0.0' )

    parser.add_option("-o", "--output", dest="outFile", default=None,
              help="Output NetCDF file (default: same was input with .cam. replaced with .camp.)")

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
        cam_File = args[0]
    else:
        parser.error("must have 1 argument: cam_File")

    if options.Vars is not None:
        options.Vars = options.Vars.split(',')

    if options.outFile is None:
        options.outFile = cam_File.replace('.cam.','.camp.')
    if options.outFile == cam_File:
        options.outFile = 'camp.nc4' # do not overwite input
        
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
    os.system("cp "+cam_File+" interm.nc")

    # Open the input file
    # -------------------
    nc0 = Dataset(cam_File, "r", format=options.format)
    nc = Dataset("interm.nc", 'a', format=options.format)
    # Air density  
    # -----------
    append_density(nc, zlib=False)

    # Time range
    # ----------
    tyme = getTyme(nc0)
    print(tyme)


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
 
    # append other vars 
    append_vars(nc, options, zlib=False)
