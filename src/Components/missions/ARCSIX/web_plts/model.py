# Module of functions to manipulate GEOS/NNR MODIS Level 3 files
import numpy as np
import numpy.ma as ma
import os.path
from netCDF4 import Dataset
import calendar
import copy
import datetime

# Need trap that if not found return NaN
def read1(filename,varn='SUEXTCOEF',z=-1):
    lonname = 'lon'
    latname = 'lat'
    levname = 'lev'

    if os.path.isfile(filename)==True:
       rc = 0
       with Dataset(filename,'r') as ncid:
           lon = ncid.variables[lonname][:]
           lat = ncid.variables[latname][:]
           if z < 1:
               try:
                   ncid.variables[levname][:]
                   lev = ncid.variables[levname][:]
                   ext = np.squeeze(ncid.variables[varn][:])
               except:
                   lev = 1
                   ext = np.squeeze(ncid.variables[varn][:])
           else:
               lev = ncid.variables[levname][:][z]
               ext = np.squeeze(ncid.variables[varn][:][:,z,:,:])
    else:
        rc = 1
        ext = 1.
        lon = 1.
        lat = 1.
        lev = 1.

    return ext, lon, lat, lev, rc

# Read a zonal average
def readzonal(filename,varn='SUEXTCOEF',z=-1):
    print(varn,expid)
    ext,lon,lat,levs,rc = readday(filename,z=z,varn=varn)
    print(ext.shape,rc,lev)
    if rc == 0:
        print(ext.shape)
        if z < 0:
            extz = ma.mean(ext,axis=1)
        else:
            extz = ma.mean(ext,axis=2)
    else:
        extz = 1
        lat = 1
    return extz, lat, rc
