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
import collections

from curv2llp   import *

def get_cP(n,ncin):
    """
    Given time level "n", return 3D pressure for this time.
    """
    lon = ncin.variables['lon']
    lat = ncin.variables['lat']
    a = ncin.variables['hyam'][:]
    b = ncin.variables['hybm'][:]
    nx, ny, nz = lon.shape[0], lat.shape[0], len(a)
    ps = ncin.variables['PS'][n,:,:]
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
        lon = ncin.variables['lon']
        lat = ncin.variables['lat']
        a = ncin.variables['hyam'][:]
        b = ncin.variables['hybm'][:]
        nx, ny, nz = lon.shape[0], lat.shape[0], len(a)
        ps = ncin.variables['PS'][n,:,:]
        cP = npy.zeros((nz,ny,nx))
        for k in range(nz):
            cP[k,:,:] = ( a[k] + b[k]*ps[:,:] )/100.0   # hPa
        return cP

def get_cPi(n,ncin):
    """
    Given time level "n", return 3D pressure (at model interfaces) for this time.
    """
    lon = ncin.variables['lon']
    lat = ncin.variables['lat']
    a = ncin.variables['hyai'][:]
    b = ncin.variables['hybi'][:]
    nx, ny, nz = lon.shape[0], lat.shape[0], len(a)
    ps = ncin.variables['PS'][n,:,:]
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

    tk = ncin.variables['T'][:,:,:,:]
    rho = tk
    nt, nz, ny, nx = tk.shape
    
    for n in range(nt):
        print n
        cP = get_cP(n,ncin)
        Rd = 287.05 
        for k in range(nz):
            rho[n,k,:,:] = cP[k,:,:]*100.0/Rd/tk[n,k,:,:]

    return rho
#---
def append_lml(ncin, options, rho, zlib=False):
    """
    Append vars at lowest model level to netcdf file.
    """

    # Open output NC file
    # -------------------
    nc = Dataset(options.outFile,'a',format=options.format)

    for v in ['NH4','so4_a1','so4_a2','so4_a3','pom_a4',
       'soa_a1','soa_a2','bc_a1','bc_a4','dst_a1',
       'dst_a2','dst_a3','ncl_a1','ncl_a2','ncl_a3','CO','NO2','NO','O3','SO2']:
        this = nc.createVariable(v+'_lml','f4',('time','lat','lon'),zlib=zlib)
        var = ncin.variables[v][:,:,:,:]
        this.long_name = ''
        this.missing_value = UNDEF # out of domain
        this.units = ''
        nt, nz, ny, nx = var.shape
        for n in range(nt):
            this[n,:,:] = var[n,87,:,:]

    v = 'rho'
    this = nc.createVariable(v+'_lml','f4',('time','lat','lon'),zlib=zlib)
    this.long_name = ''
    this.missing_value = UNDEF # out of domain
    this.units = ''
    nt, nz, ny, nx = rho.shape
    for n in range(nt):
        this[n,:,:] = rho[n,87,:,:]

#   v = 'rho'
#   this = nc.createVariable(v,'f4',('time','lev','lat','lon'),zlib=zlib)
#   this.long_name = ''
#   this.missing_value = UNDEF # out of domain
#   this.units = ''
#   nt, nz, ny, nx = rho.shape
#   print this.shape, rho.shape
#   for n in range(nt):
#       this[n,:,:,:] = rho[n,:,:,:]


    nc.close()

#---
def append_vars(ncin, options, rho, zlib=False):
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

    for v in ['pom_a4','soa_a1','soa_a2']:
        oc = oc + ncin.variables[v][:,:,:,:]/1.4 # 1.4 is the ratio of pom emissions vs oc emissions

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

    # total column pm2.5
    # ---------------
    pm25 = ncin.variables['PM25'][:,:,:,:]  # kg m-3
    this = nc.createVariable('tcpm25','f4',('time', 'lat', 'lon'),zlib=zlib)
    this.long_name = 'total column PM2.5'
    this.missing_value = UNDEF # out of domain
    this.units = 'kg m-2'
    for n in range(nt): 
        tc = npy.sum(pm25[n,:,:,:]*delp[n,:,:,:]*100.0/rho[n,:,:,:]/g, axis=0)
        this[n,:,:] = tc[:,:]

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
    rP = '1000,975,950,925,900,875,850,825,800,775,750,725,700,650,600,550,500,450,400,350,300,250,200,150,100'

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


    # Open the input file
    # -------------------
    nc0 = Dataset(cam_File, "r", format=options.format)

    # Air density  
    # -----------
    rho = append_density(nc0, zlib=False)

    # Time range
    # ----------
    tyme = getTyme(nc0)
    print(tyme)


    # Coordinates
    # -----------
    lon = nc0.variables['lon'][:]
    lat = nc0.variables['lat'][:]

    # Instantiate regridding class
    # ----------------------------
    print 'Start regridding'
    r = myCurv2LLP(lon,lat,options.rP)
    print 'FInished regridding'
    
    # Write output file
    # -----------------
    r.writeNC ( tyme, options, nc0, zlib=False )
    
    # append variables at the lowest model level 
    append_lml(nc0, options, rho, zlib=False)
 
    # append other vars 
    append_vars(nc0, options, rho, zlib=False)
