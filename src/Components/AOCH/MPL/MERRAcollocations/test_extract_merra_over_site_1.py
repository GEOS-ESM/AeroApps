#!/usr/bin/env python3
"""
Mar/21/2025 
"""

from datetime import datetime
import netCDF4 as nc

import numpy as np
import sys, os

from pyobs.sampler import TRAJECTORY   #for extracting MERRA along track
from pyobs.sampler import STATION   #for extracting MERRA along track

from pyobs.sampler import addVertCoord #for getting MERRA-2 vertical coordinate
from pyobs.aop import G2GAOP #For getting MERRA2 AOPs 


import matplotlib
import matplotlib.pyplot as plt

import shutil ### use to copy input file and then add new arrays


def parse_datetime(date_string):
    return datetime.strptime(date_string, '%Y-%m-%d %H:%M:%S.%f')

def process_file(file_path):
    # Read the entire file
    with open(file_path, 'r') as file:
        # Read all lines, skip the header
        lines = file.readlines()[1:]

    # Initialize lists for each column
    datetimes = []
    orbit_numbers = []
    source_filenames = []
    lines_col = []
    columns_col = []

    # Process each line
    for line in lines:
        # Split the line by comma
        parts = line.strip().split(',')
        
        # Append each part to its respective list
        datetimes.append(parse_datetime(parts[0]))
        orbit_numbers.append(int(parts[1]))
        source_filenames.append(parts[2])
        lines_col.append(int(parts[3]))
        columns_col.append(int(parts[4]))

    return {
        'datetimes': datetimes,
        'orbit_numbers': orbit_numbers,
        'source_filenames': source_filenames,
        'lines': lines_col,
        'columns': columns_col
    }
    
def print_data_summary(data):
    print("First 5 datetimes:", data['datetimes'][:5])
    print("First 5 orbit numbers:", data['orbit_numbers'][:5])
    print("First 5 source filenames:", data['source_filenames'][:5])
    print("First 5 lines:", data['lines'][:5])
    print("First 5 columns:", data['columns'][:5])
    print(f"Total rows read: {len(data['datetimes'])}")

def print_specific_row(data, row_index):
    print(f"\nData for row {row_index + 1}:")
    print(f"Datetime: {data['datetimes'][row_index]}")
    print(f"Orbit number: {data['orbit_numbers'][row_index]}")
    print(f"Source filename: {data['source_filenames'][row_index]}")
    print(f"Line: {data['lines'][row_index]}")
    print(f"Column: {data['columns'][row_index]}")
#-----------------------------------------

##### Set path of output file in NoNBackup so I do not fill my Discover quota
#pth_in_sat='/discover/nobackup/sgasso/Files/Satellite/SliceTROCAL/' # INPUT w/CALIOP-TROP collocations
#pth_in_mod='/discover/nobackup/sgasso/Files/Satellite/MERRA2/'      # INPUT w/MERRA2-CAL collocations
#pth_out   ='/discover/nobackup/sgasso/Files/Satellite/Triple/'      # OUPUT w/MERRA-CAL-TROP collocations

# Name_MPL_sites , name,lat, lon, altitude(km), time1,time2
# GSFC                ,38.9930,-76.8400,0.050 , 17:05,19:05
# Santa_Cruz_Tenerife ,28.4720,-16.2470,0.052  ,13:15,15:15

############# Set input file with CAL&TROP collocation
file_path  ='test_input_trop_overpass_1.txt'# 
site_xy = [38.9930,-76.8400]

### Output name for MERRA2 data along CAL track

# out_MERRA = 'm2_calipso_sampled_'+CAL_IOWA_File[12:-3]+'.nc'
out_MERRA = 'm2_test.nc'
save_merra='False'

#print('Hello')
### Now set output filename for triple collocation 
#output_filename = 'AlongTrack_MERRA-IOWA-CALIOP_' + caldata.cal_filename[35:-15] + '_3.1.nc'
#output_filename = 'AlongTrack_MERRA-IOWA-CALIOP_' + CAL_IOWA_File[12:-3] + '.nc'


# Process the file
data = process_file(file_path)

## Print a summary of the data
# print_data_summary(data)
## Print information for a specific row (e.g., the third row)
# print_specific_row(data, 2)
## You can now use the data in your script
## For example, to get all datetimes:
# all_datetimes = data['datetimes']
## Or to get the orbit number of the 5th row:
# fifth_orbit_number = data['orbit_numbers'][4]

#sys.exit()

####
#### Process CALIPSO data and compute weighted height
####
npts=1
# npts = len(caldata.cal_lat) # number of profiles
# nlev = len(caldata.cal_layer_km) # number of layer
# # Get time coordinate
# time  = caldata.time
# ntime = len(time)
# time  = np.repeat(time.reshape(ntime,1),nlev,axis=1)

###### CALIOP ext and z : set negative values to zero so they do not contriute to the column integration
### # have to flip order of array so it agrees with bext, for some reason both bext and z are not in the same order.
# cal_z = caldata.cal_layer_km[::-1] 
# cal_maxz=cal_z[0] # save this for comparing with MERRA
# # ilast=nlev-10+np.where(cal_z[-10:] < 0)[0][0]-1### Now find where the first non-zero height, I may not need this but it maybe useful
# # ifirst=np.where(cal_z[0:14] > 0)[0][0] # if using cal_z = caldata.cal_layer_km
# # Create intermediate array and set to 0 all negative values
# cal_ext532  = caldata.cal_ext; cal_ext532  = np.where(caldata.cal_ext < 0, 0, caldata.cal_ext) # set to 0 all negative values
# ### Compute Weighted CAL height at two wavelenghts
# # A1 =  np.sum(cal_ext532[:,0:ilast]*cal_z[0:ilast],axis=1)
# # B1 =  np.sum(cal_ext532[:,0:ilast],axis=1)
# A1 =  np.sum(cal_ext532*cal_z,axis=1)
# B1 =  np.sum(cal_ext532,axis=1)
# cal_hw1=np.divide(A1.data, B1, out=np.full_like(A1, np.nan), where=B1!=0)

#sys.exit()
# print(cal_hw[~np.isnan(cal_hw)])

#####
##### now do the same for MERRA-2
#####

### # set up some file names
m2data = 'inst3_3d_aer_Nv' # do not change this 
#m2data = np.array(['inst3_3d_aer_Nv','inst3_3d_aer_Nv'])

##### Read the CALIPSO file, get the lat, lon, times from the file
# caldata = CALIPSO_L2(calipsoFile,Verbose=True)
##### create a trajectory object
# time, lon, lat = zip(*[[data['datetimes'][0], site_xy[1], site_xy[0]],
                       # [data['datetimes'][0], site_xy[1], site_xy[0]]])
#### -----
data_list = [[data['datetimes'][0], site_xy[1], site_xy[0]],
             [data['datetimes'][1], site_xy[1], site_xy[0]]]
time_list = [item[0] for item in data_list]
lon_list  = [item[1] for item in data_list]
lat_list  = [item[2] for item in data_list]
time_array = np.array(time_list)
lon_array = np.array(lon_list)
lat_array = np.array(lat_list)

traj = TRAJECTORY(time_array,lon_array,lat_array,m2data) # sample the MERRA-2 dataset along the trajectory, and return an xarray dataset
traj2 = STATION(['GSFC','GSFC'],lon_array,lat_array,m2data,time_array)
# sample the MERRA-2 dataset along the trajectory, and return an xarray dataset
traj_ds = traj.sample()
traj2_ds = traj2.sample()
# write sampled data to a netcdf file
print("DONE extracting MERRA data ")
sys.exit()
if save_merra == 'True':
   print('\n Saving MERRA data file.... ', out_MERRA)
   # print('      In folder ', pth_in_mod)
   traj_ds.to_netcdf(out_MERRA)
   #traj_ds.to_netcdf(pth_in_mod+out_MERRA)
   print(" Done Savings!\n")


# sys.exit()


# ####### -------------------
# # read the sampled aerosol profile data and optical tables
# # set up some filenames
config = 'm2_pm25.yaml' # this configuration file can be found in src/config
# MERRA2_File = 'm2_calipso_sampled.nc4'
optics = G2GAOP(out_MERRA,config=config)

# Get BEXT at 532 nm and 1064nm
# breakpoint()
ext532  = optics.getAOPext(wavelength=532)
#ext1064 = optics.getAOPext(wavelength=1064)

# # Get MERRA-2 vertical coordinate
ext532.pipe(addVertCoord)

# # Get time coordinate
time  = ext532.time.values
ntime = ext532.dims['time']
nlev  = ext532.dims['lev']
time  = np.repeat(time.reshape(ntime,1),nlev,axis=1)

# ## Sanity check. Just to make sure the ingested EXT are same lenght as CALIPSO
# if npts != ext532.EXT.shape[0] : sys.exit('Number of points from MERRA do not match Npixels from CALIPSO')

# ### Now compute MERRA weighted height
me_zlayer =ext532.Z[0].values   # in km

# # Find the index of the closest element to Maximuing CALIOP height
# indx = np.argmin(np.abs(me_zlayer - cal_maxz))
# me_zlayer2=np.tile(me_zlayer[indx:], (npts, 1))
# me_text532 =ext532.EXT.values[:,indx:] # in km-1
# #me_text1064=ext1064.EXT.values[:,indx:] # in km-1

# Find the index of the closest element to Maximuing CALIOP height
me_zlayer2=np.tile(me_zlayer, (npts, 1))
me_text532 =ext532.EXT.values # in km-1
#me_text1064=ext1064.EXT.values[:,indx:] # in km-1


if len(np.where(me_text532<0)[0])>1 : sys.exit('Some EXT are <0 in MERRA. Check')
A3=np.sum(me_zlayer2 * me_text532,axis=1)
B3=np.sum(me_text532,axis=1)
me_hw532=np.divide(A3 , B3) # km
#A4=np.sum(me_zlayer2 * me_text1064,axis=1)
#B4=np.sum(me_text1064,axis=1)
#me_hw1064=np.divide(A4 , B4) # km

# ##### Done computing weighted heights

sys.exit()

# #### NOW save output
# # Create the new file directly in the output directory
# # output_filename = 'AlongTrack_MERRA-IOWA-CALIOP_' + caldata.cal_filename[35:-15] + '.nc'
# output_file = os.path.join(pth_out, output_filename)  # This creates a proper path string
# # Copy the original file directly to the output location
# shutil.copy2(os.path.join(pth_in_sat, CAL_IOWA_File), output_file)
# ### Put together array to add to the output file
# #combined_array = np.array([cal_hw1,cal_hw2, caldata.tro_height, me_hw532,me_hw1064]).transpose()
# combined_array = np.array([cal_hw1,caldata.tro_height, me_hw532]).transpose()
# # sys.exit()
# # Now open the file in the output directory
# with nc.Dataset(output_file, 'a') as dataset:
    # # Create a new dimension for the number of valid data points
    # if 'nlines' not in dataset.dimensions:
        # dataset.createDimension('nlines', combined_array.shape[0])
    # # Create a new dimension for the height types
    # if 'height_types' not in dataset.dimensions:
        # dataset.createDimension('height_types', combined_array.shape[1])
    
    # # Create the new variable
    # combined_var = dataset.createVariable('combined_heights', np.float32, ('nlines', 'height_types'))
    
    # # Add attributes to describe the variable
    # combined_var.long_name = 'Aerosol Heights from IOWA, MERRA and CALIOP along CALIOP Track'
    # combined_var.units = 'km'
    # # combined_var.description = 'Column 0,1: CALIPSO weighted height 532 and 1064nm, Column 2: IOWA AOCH , Column 3-4: MERRA-2 weighted height 531,1064'
    # combined_var.description = 'Column 0,1: CALIPSO weighted height 532 , Column 2: IOWA AOCH , Column 3-4: MERRA-2 weighted height 531,1064'
    
    # # Write the data to the variable
    # combined_var[:] = combined_array
# print(f"Combined height data has been added to \n {pth_out}{output_filename}")
