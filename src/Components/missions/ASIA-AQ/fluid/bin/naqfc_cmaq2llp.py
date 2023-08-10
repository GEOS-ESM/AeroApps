#!/usr/bin/env python
#
# Convert NAQFC-CMAQ output to COARDS compliant file on a regular lat-lon-pressure grid.
#

from dateutil.parser import parse as isoparser
from optparse   import OptionParser
from datetime   import datetime, timedelta
from grads      import GrADS
from netCDF4    import Dataset, MFDataset
import numpy as npy
import os

from curv2llp   import *

#---
def getTyme(ncin):
    t = ncin.variables['TFLAG'][:,0,:] #TFLAG:units = "<YYYYDDD,HHMMSS>" ;
    nt, ns = t.shape
    day = timedelta(days=1)
    tyme = npy.zeros(nt,dtype=datetime)

    for i in range(nt):
        dstr = str(t[i,0]).zfill(7)+str(t[i,1]).zfill(6)
        tyme[i] = datetime(int(dstr[0:4]),1,1,int(dstr[7:9])) + (int(dstr[4:7])-1)*day
    return tyme

#---
class myCurv2LLP(Curv2LLP):
    
    def get_cP(self,n,ga):
        """
        Given time level "n", return 3D pressure for this time.
        """
        ga('set t %d'%n)
        fh = ga.query('file')
        nz, ny, nx = fh.nz, fh.ny, fh.nx
        ga('set z 1 %d'%fh.nz)
        ga('set y 1 %d'%fh.ny)
        ga('set x 1 %d'%fh.nx)
        p = ga.expr('pres')
        p = p/100.0   # in hPa
        print(p.shape)
        return p 
    #def get_cP(self,n,ncin):
    #    return ncin.variables['PRES'][n,:,:,:]/100.  # in hPa

#---
def getCoords(ga):

    fh = ga.query('file')

    ga('set x 1 %d'%fh.nx)
    ga('set y 1')
    lon = ga.expr('lon')
    
    ga('set x 1')
    ga('set y 1 %d'%fh.ny)
    lat = ga.expr('lat')

    return (lat,lon)

#---
def append_vars(ncin, options, zlib=False):
    """
    Append total column vars to output file 
    """
    # Open output NC file
    # -------------------
    nc = Dataset(options.outFile,'a',format=options.format)

    p = ncin.variables['PRES'][0,:,:,:]  # Pa
    nz, ny, nx = p.shape 
    oc1 = ncin.variables['AOLGAJ'][0,:,:,:]   # ug m-3 
    oc2 = ncin.variables['AOLGBJ'][0,:,:,:]
    oc3 = ncin.variables['AORGCJ'][0,:,:,:]
    co = ncin.variables['CO'][0,:,:,:]*1e-6*28.0/28.97 # ppm -> kg kg-1
    z = ncin.variables['ZF'][0,:,:,:]  # full layer height above ground in m 
    dz = npy.zeros((nz,ny,nx))
   
    for n in range(nz):
        if n==0:
            dz[n,:,:] = z[n,:,:]
        else:
            dz[n,:,:] = z[n,:,:] - z[n-1,:,:]
    
    tk = ncin.variables['TA'][0,:,:,:]
    rho = npy.zeros((nz,ny,nx))
    Rd = 287.05 
    rho[:,:,:] = p[:,:,:]/Rd/tk[:,:,:]
    

    # total column OC 
    # ---------------
    this = nc.createVariable('tcoc','f4',('time', 'lat', 'lon'),zlib=zlib)
    this.long_name ='total cloumn OC'
    this.missing_value = UNDEF # out of domain
    this.units = 'kg m-2'
    this[0,:,:] =npy.sum((oc1 + oc2 + oc3)*1e-9*dz, axis=0) 

    # total column CO
    # ---------------
    this = nc.createVariable('tcco','f4',('time', 'lat', 'lon'),zlib=zlib)
    this.long_name ='total cloumn CO'
    this.missing_value = UNDEF # out of domain
    this.units = 'kg m-2'
    this[0,:,:] =npy.sum( co*rho*dz, axis=0) 


#---        
if __name__ == "__main__":

    format = 'NETCDF4'
    algo = 'linear'
    rP = '925,850,700,600'

#   Parse command line options
#   --------------------------
    parser = OptionParser(usage="Usage: %prog [OPTIONS] cmaq_file",
                          version='1.0.0' )

    parser.add_option("-o", "--output", dest="outFile", default=None,
              help="Output NetCDF file (default: same as input, end with .llp.nc4)")

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
        cmaq_file = args[0]
    else:
        parser.error("must have 1 argument: cmaq_file")

    if options.Vars is not None:
        options.Vars = options.Vars.split(',')

    if options.outFile is None:
        dirn = os.path.dirname(cmaq_file)
        base = os.path.basename(cmaq_file)
        options.outFile = dirn + '/' + base.replace('.ncf','.llp.nc4')
    if options.outFile == cmaq_file:
        options.outFile = 'cmaq_llp.nc4' # do not overwite input
        
    options.rP = [ float(p) for p in options.rP.split(',') ]
    
    if options.ignore is not None:
        options.ignore = options.ignore.split(',')
    else:
        options.ignore = []
        
    options.ignore += ['pres'] # these are coordinates

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
    nc = Dataset(cmaq_file)

    # Time range
    # ----------
    tyme = getTyme(nc)

    # Get lat lon 
    # ----------
    ga = GrADS(Window=False,Echo=False)
    cmaq_ctl = cmaq_file.replace('.nc', '.ctl')
    fh = ga.open(cmaq_ctl)
    lat, lon = getCoords(ga)

    # Instantiate regridding class
    # ----------------------------
    r = myCurv2LLP(lon,lat,options.rP)

    # Write output file
    # -----------------
    r.writeNC ( tyme, options, ga, zlib=False )
    
    # Append vars 
    # -----------
    append_vars(nc, options, zlib=False)


    
