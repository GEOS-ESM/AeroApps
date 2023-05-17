#!/usr/bin/env python3
"""
Makes netcdf files of input data for NN trainig dataset
"""
import os
from dateutil.parser import parse
from netCDF4 import Dataset
import numpy as np
from math import cos, radians
#---
if __name__ == '__main__':
  outdir = 'vnncLUT'
  filename = 'LUT_angles_wind.nc4'


 
  sza     = np.linspace(0.15,1,18)
  vza     = np.linspace(0.15,1,18)
  sza     = np.degrees(np.arccos(sza))
  vza     = np.degrees(np.arccos(vza))
  raa     = np.degrees(np.arccos([-1, -0.5, 0, 0.5, 1]))  

  u10m    = np.linspace(0,20,5)
  v10m    = np.array([0])
  u10m[0] = 0.10
  v10m[0] = 0.10

  # create file
  outfile        = outdir + '/' + filename
  ncOut          = Dataset(outfile,'w',format='NETCDF4_CLASSIC')
  ncOut.comment  = 'angles and wind speeds for ocean NN training data'

  # create dimensions
  sza_dim      = ncOut.createDimension('sza',len(sza))
  vza_dim      = ncOut.createDimension('vza',len(vza))
  raa_dim      = ncOut.createDimension('raa',len(raa))  
  uwind_dim     = ncOut.createDimension('uwind',len(u10m))  
  vwind_dim     = ncOut.createDimension('vwind',len(v10m))    

  # write out variables
  varobj = ncOut.createVariable('sza','f4',('sza',))
  varobj.long_name = 'solar zenith angle'
  varobj.units     = 'degrees'
  varobj[:] = sza

  varobj = ncOut.createVariable('vza','f4',('vza',))
  varobj.long_name = 'viewing zenith angle'
  varobj.units     = 'degrees'
  varobj[:] = vza 

  varobj = ncOut.createVariable('raa','f4',('raa',))
  varobj.long_name = 'relative azimuth angle'
  varobj.units     = 'degrees'
  varobj[:] = raa   

  varobj = ncOut.createVariable('u10m','f4',('uwind',))
  varobj.long_name = 'u10m'
  varobj.units     = 'm/2'
  varobj[:] = u10m

  varobj = ncOut.createVariable('v10m','f4',('vwind',))
  varobj.long_name = 'v10m'
  varobj.units     = 'm/2'
  varobj[:] = v10m

  ncOut.close()
