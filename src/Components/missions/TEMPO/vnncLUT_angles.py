#!/usr/bin/env python
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
  filename = 'LUT_angles_albedo.nc4'

  sza_space    = 0, 80, 9
  vza_space    = 0, 80, 9 
  raa_space    = 0, 180, 5
  albedo_space1 = 0, 0.1, 5
  albedo_space2 = 0.2, 1, 5

  # sza_space    = 5, 80, 2
  # vza_space    = 5, 80, 2
  # raa_space    = 30, 30, 1
  # albedo_space1 = 0.1, 0.1, 1 
  # albedo_space2 = 1, 1, 1

  sza     = np.linspace(cos(radians(sza_space[0])),cos(radians(sza_space[1])),sza_space[2])
  vza     = np.linspace(cos(radians(vza_space[0])),cos(radians(vza_space[1])),vza_space[2])
  # raa     = np.linspace(raa_space[0],raa_space[1],raa_space[2])

  sza     = np.array([0.15, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.0])
  vza     = np.array([0.15, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.0])
  # raa     = np.degrees(np.arccos([-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1]))
  raa     = np.degrees(np.arccos([-1, -0.5, 0, 0.5, 1]))  

  sza     = np.array([0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.60, 0.70, 0.80, 0.90, 0.95, 1.0])
  vza     = np.array([0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.60, 0.70, 0.80, 0.90, 0.95, 1.0])
  # raa     = np.degrees(np.arccos([-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1]))
  raa     = np.degrees(np.arccos([-1, -0.5, 0, 0.5, 1]))  

  sza     = np.linspace(0.15,1,18)
  vza     = np.linspace(0.15,1,18)

  albedo  = np.append(np.linspace(albedo_space1[0],albedo_space1[1],albedo_space1[2]),np.linspace(albedo_space2[0],albedo_space2[1],albedo_space2[2]))
  nalbedo = albedo_space1[2] + albedo_space2[2]

  sza     = np.degrees(np.arccos(sza))
  vza     = np.degrees(np.arccos(vza))

  # create file
  outfile        = outdir + '/' + filename
  ncOut          = Dataset(outfile,'w',format='NETCDF4_CLASSIC')
  ncOut.comment  = 'angles and albedo for NN training data'

  # create dimensions
  sza_dim      = ncOut.createDimension('sza',len(sza))
  vza_dim      = ncOut.createDimension('vza',len(vza))
  raa_dim      = ncOut.createDimension('raa',len(raa))  
  albedo_dim   = ncOut.createDimension('albedo',len(albedo))    

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

  varobj = ncOut.createVariable('albedo','f4',('albedo',))
  varobj.long_name = 'surface reflectance'
  varobj.units     = 'none'
  varobj[:] = albedo

  ncOut.close()
