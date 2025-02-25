#!/usr/bin/env python3
"""
Jan/2/2025 - _ver _3.2, it is base in 3.1 but adapted to run a batch of files as opposed to 1 in
ver 3.1 
Dec/30/2024 - ver _3.1 , this code just runs that cases after AGU. 

Dec/30/2024 - _3 , this version works Discover and reproduces cases used in AGU.
A couple of things, 
1) unlike Calculon version, it incorporated the code extract_merra_along_caliop_2.1.py
so now I only need to run CALTROP_overlap_4.1py and this code.
2) I had to copy the calipso_l2.py and aop.py codes from Calculon , the standard version fron 
Discover did not have a couple of operations needed for extracting data from MERRA.

Dec/16/2024 - same as collocated_MERRA-CALIOP-IOWA_2.py in my Python/Discover folder
but it includes the script /home/sgasso/Python/AOCH/extract_merra_along_caliop_1.py from discover

NEED to finish it!

Dec/2/2024 - Now works and running it for cases for AGU

Inputs required 
CAL_IOWA_File which is the output of CALTROP_overlap_3.py 
(currently running in Calculon)
MERRA2_File : output of extract_merra_along_caliop_1.py, eventually this should be merged with this code.


"""

from datetime import datetime
import netCDF4 as nc

import numpy as np
import sys, os

from pyobs.sampler import TRAJECTORY   #for extracting MERRA along track
from pyobs.sampler import addVertCoord #for getting MERRA-2 vertical coordinate
from pyobs.aop import G2GAOP #For getting MERRA2 AOPs 


import matplotlib
import matplotlib.pyplot as plt

import shutil ### use to copy input file and then add new arrays


class CALDATA:
    def __init__(self, **kwargs):
        for key, value in kwargs.items():
            setattr(self, key, value)

def load_caldata(file_path, variables):
    caldata = CALDATA()
    with nc.Dataset(file_path, 'r') as dataset:
        for var in variables:
            if var in dataset.variables:
                if var == 'cal_time':
                    # Special handling for cal_time
                    time_var = dataset.variables[var]
                    times = time_var[:]
                    time_units = time_var.units
                    time_calendar = time_var.calendar
                    datetime_objects = nc.num2date(times, units=time_units, calendar=time_calendar)
                    setattr(caldata, 'time', np.array([datetime(dt.year, dt.month, dt.day, 
                                                dt.hour, dt.minute, dt.second) 
                                       for dt in datetime_objects], dtype=object))
                elif var == 'calipso_filename' :
                    # Special handling for calipso_filename
                    setattr(caldata, 'cal_filename', dataset.variables[var][0])
                elif var == 'iowa_filename' :
                    # Special handling for IOWA_filename
                     setattr(caldata, 'iowa_filename', dataset.variables[var][0])
                else:
                    # General case for other variables
                    data = dataset.variables[var][:]
                    if data.ndim == 1:
                        setattr(caldata, var, data[:])
                    elif data.ndim == 2:
                        setattr(caldata, var, data[:, :])
                    else:
                        print(f"Warning: Variable '{var}' has unexpected dimensions. Skipping.")
            else:
                print(f"Warning: Variable '{var}' not found in the NetCDF file.")
 
    return caldata

#-----------------------------------------

##### Set path of output file in NoNBackup so I do not fill my Discover quota
pth_in_sat='/discover/nobackup/sgasso/Files/Satellite/SliceTROCAL/' # INPUT w/CALIOP-TROP collocations
pth_in_mod='/discover/nobackup/sgasso/Files/Satellite/MERRA2/'      # INPUT w/MERRA2-CAL collocations
pth_out   ='/discover/nobackup/sgasso/Files/Satellite/Triple/'      # OUPUT w/MERRA-CAL-TROP collocations

############# Set input file with CAL&TROP collocation
CAL_IOWA_File ='IOWA-CALIOP_2020-09-10T21-04-23ZD_AGU.nc'# smoke  case
#CAL_IOWA_File ='IOWA-CALIOP_2021-08-28T21-07-48ZD_case2.nc' # 
#CAL_IOWA_File ='IOWA-CALIOP_2021-08-28T14-33-47ZD_AGU.nc' # dust case


var_list = ['calipso_filename', 'cal_lat', 'cal_lon', 'cal_ext', 'cal_ext2','cal_tback', 'cal_time','cal_layer_km',\
            'tro_lat','tro_lon','tro_height','iowa_filename','tro_height','tro_aod680']
print('Loading data from : ',pth_in_sat + CAL_IOWA_File)
caldata = load_caldata(pth_in_sat + CAL_IOWA_File, var_list)

### Output name for MERRA2 data along CAL track
out_MERRA = 'm2_calipso_sampled_'+CAL_IOWA_File[12:-3]+'.nc'
### Now set output filename for triple collocation 
#output_filename = 'AlongTrack_MERRA-IOWA-CALIOP_' + caldata.cal_filename[35:-15] + '_3.1.nc'
output_filename = 'AlongTrack_MERRA-IOWA-CALIOP_' + CAL_IOWA_File[12:-3] + '.nc'

### Dec/30/ note: saving MERRA file is taking a lot of time and I need it only once per case
### So, first check if the file already exists, if no, go ahead read and save.
### if yes, do not save.
### 
if os.path.isfile(pth_in_mod+out_MERRA):
   print("\n There is a file with same name, so MERRA data will be extracted but not saved\n")
   save_merra='False'
else:
   print("\nNo MERRA file found, so file will be saved. It will take about 1 min. \n")
   save_merra='True'

#### Now with the CAL lat , lon go to MERRA and extract the slice of data
#### list so basic info
# for attr_name in dir(caldata):
  # if not attr_name.startswith('__'):  # Skip built-in attributes
      # attr_value = getattr(caldata, attr_name)
      # if isinstance(attr_value, np.ndarray):
          # print(f"Variable: {attr_name}")
          # print(f"  Dimensions: {attr_value.shape}")
          # print(f"  Data type: {attr_value.dtype}")
      # elif attr_name.endswith('filename'):
          # print(f"Variable: {attr_name}")
          # print(f"  Value: {attr_value}")
      # print("--------------------")

####
#### Process CALIPSO data and compute weighted height
####
npts = len(caldata.cal_lat) # number of profiles
nlev = len(caldata.cal_layer_km) # number of layer
# Get time coordinate
time  = caldata.time
ntime = len(time)
time  = np.repeat(time.reshape(ntime,1),nlev,axis=1)

###### CALIOP ext and z : set negative values to zero so they do not contriute to the column integration
### # have to flip order of array so it agrees with bext, for some reason both bext and z are not in the same order.
cal_z = caldata.cal_layer_km[::-1] 
cal_maxz=cal_z[0] # save this for comparing with MERRA
# ilast=nlev-10+np.where(cal_z[-10:] < 0)[0][0]-1### Now find where the first non-zero height, I may not need this but it maybe useful
# ifirst=np.where(cal_z[0:14] > 0)[0][0] # if using cal_z = caldata.cal_layer_km
# Create intermediate array and set to 0 all negative values
cal_ext532  = caldata.cal_ext; cal_ext532  = np.where(caldata.cal_ext < 0, 0, caldata.cal_ext) # set to 0 all negative values
cal_ext1064 = caldata.cal_ext2;cal_ext1064 = np.where(caldata.cal_ext2 < 0, 0, caldata.cal_ext2) # set to 0 all negative values
### Compute Weighted CAL height at two wavelenghts
# A1 =  np.sum(cal_ext532[:,0:ilast]*cal_z[0:ilast],axis=1)
# B1 =  np.sum(cal_ext532[:,0:ilast],axis=1)
A1 =  np.sum(cal_ext532*cal_z,axis=1)
B1 =  np.sum(cal_ext532,axis=1)
cal_hw1=np.divide(A1.data, B1, out=np.full_like(A1, np.nan), where=B1!=0)

# A2 =  np.sum(cal_ext1064[:,0:ilast]*cal_z[0:ilast],axis=1)
# B2 =  np.sum(cal_ext1064[:,0:ilast],axis=1)
A2 =  np.sum(cal_ext1064*cal_z,axis=1)
B2 =  np.sum(cal_ext1064,axis=1)
cal_hw2=np.divide(A2.data, B2, out=np.full_like(A2, np.nan), where=B2!=0)


# sys.exit()
# print(cal_hw[~np.isnan(cal_hw)])

#####
##### now do the same for MERRA-2
#####

#### First sample the model along the track and save it in intermediate file (for now)
### now from file extract_merra_along_caliop_21.py
### Variables to read
# var_cal = ['calipso_filename', 'cal_lat', 'cal_lon','cal_time']
# print('   Reading MERRA file',pth_in_sat+caltro_File)
# caltro_data = load_caldata(pth_in_sat+caltro_File,var_cal)
# caltro_data = load_caldata(pth_in_sat+caltro_File)
##### Misc print infor vars used in model sampling.
# print("CALIPSO Filename:", caldata.cal_filename)
# print("Latitude   shape:", caldata.cal_lat.shape)
# print("Longitude  shape:", caldata.cal_lon.shape)
# print("Time       shape:", caldata.time.shape)
### # set up some file names
m2data = 'inst3_3d_aer_Nv' # do not change this 
##### Read the CALIPSO file, get the lat, lon, times from the file
# caldata = CALIPSO_L2(calipsoFile,Verbose=True)
##### create a trajectory object
time, lon, lat = caldata.time, caldata.cal_lon[:], caldata.cal_lat[:]
traj = TRAJECTORY(time,lon,lat,m2data)
# sample the MERRA-2 dataset along the trajectory, and return an xarray dataset
traj_ds = traj.sample()
# write sampled data to a netcdf file
print("DONE extracting MERRA data extracted and collocated along CALIOP")
if save_merra == 'True':
   print('\n Saving MERRA data file.... ', out_MERRA)
   print('      In folder ', pth_in_mod)
   traj_ds.to_netcdf(pth_in_mod+out_MERRA)
   print(" Done Savings!\n")


#sys.exit()


####### -------------------
# read the sampled aerosol profile data and optical tables
# set up some filenames
config = 'm2_pm25.yaml' # this configuration file can be found in src/config
# MERRA2_File = 'm2_calipso_sampled.nc4'
optics = G2GAOP(pth_in_mod + out_MERRA,config=config)


# Get BEXT at 532 nm and 1064nm
# breakpoint()
ext532  = optics.getAOPext(wavelength=532)
ext1064 = optics.getAOPext(wavelength=1064)

# Get MERRA-2 vertical coordinate
ext532.pipe(addVertCoord)

# Get time coordinate
time  = ext532.time.values
ntime = ext532.dims['time']
nlev  = ext532.dims['lev']
time  = np.repeat(time.reshape(ntime,1),nlev,axis=1)

## Sanity check. Just to make sure the ingested EXT are same lenght as CALIPSO
if npts != ext532.EXT.shape[0] : sys.exit('Number of points from MERRA do not match Npixels from CALIPSO')

### Now compute MERRA weighted height
me_zlayer =ext532.Z[0].values   # in km
# Find the index of the closest element to Maximuing CALIOP height
indx = np.argmin(np.abs(me_zlayer - cal_maxz))
me_zlayer2=np.tile(me_zlayer[indx:], (npts, 1))
me_text532 =ext532.EXT.values[:,indx:] # in km-1
me_text1064=ext1064.EXT.values[:,indx:] # in km-1
if len(np.where(me_text532<0)[0])>1 : sys.exit('Some EXT are <0 in MERRA. Check')
A3=np.sum(me_zlayer2 * me_text532,axis=1)
B3=np.sum(me_text532,axis=1)
me_hw532=np.divide(A3 , B3) # km
A4=np.sum(me_zlayer2 * me_text1064,axis=1)
B4=np.sum(me_text1064,axis=1)
me_hw1064=np.divide(A4 , B4) # km

##### Done computing weighted heights

# sys.exit()

#### NOW save output
# Create the new file directly in the output directory
# output_filename = 'AlongTrack_MERRA-IOWA-CALIOP_' + caldata.cal_filename[35:-15] + '.nc'
output_file = os.path.join(pth_out, output_filename)  # This creates a proper path string
# Copy the original file directly to the output location
shutil.copy2(os.path.join(pth_in_sat, CAL_IOWA_File), output_file)
### Put together array to add to the output file
combined_array = np.array([cal_hw1,cal_hw2, caldata.tro_height, me_hw532,me_hw1064]).transpose()
# sys.exit()
# Now open the file in the output directory
with nc.Dataset(output_file, 'a') as dataset:
    # Create a new dimension for the number of valid data points
    if 'nlines' not in dataset.dimensions:
        dataset.createDimension('nlines', combined_array.shape[0])
    # Create a new dimension for the height types
    if 'height_types' not in dataset.dimensions:
        dataset.createDimension('height_types', combined_array.shape[1])
    
    # Create the new variable
    combined_var = dataset.createVariable('combined_heights', np.float32, ('nlines', 'height_types'))
    
    # Add attributes to describe the variable
    combined_var.long_name = 'Aerosol Heights from IOWA, MERRA and CALIOP along CALIOP Track'
    combined_var.units = 'km'
    combined_var.description = 'Column 0,1: CALIPSO weighted height 532 and 1064nm, Column 2: IOWA AOCH , Column 3-4: MERRA-2 weighted height 531,1064'
    
    # Write the data to the variable
    combined_var[:] = combined_array
print(f"Combined height data has been added to \n {pth_out}{output_filename}")
