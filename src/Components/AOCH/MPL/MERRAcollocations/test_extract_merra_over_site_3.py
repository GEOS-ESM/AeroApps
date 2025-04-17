#!/usr/bin/env python3
"""
Test code that ingests tropomi & mpl collocations, check date and location of site and pulls the MERRA data for respective
model grid cell. Then it computes AOCH  and saves in output file. 
-------------------------------------
Apr/3/2025 test_extract_merra_over_site_2.py .This code is based on ver *_1
1) it incorporates a loop to go over all dates according to input file. 
2) it save calculated MERRA AOCh in output file.
Local folder /discover/nobackup/sgasso/PROJECTS/AOCH

Mar/21/2025 
"""

from datetime import datetime, timedelta
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

#-----------------------------------------

##### Set path of output file in NoNBackup so I do not fill my Discover quota
#pth_in_sat='/discover/nobackup/sgasso/Files/Satellite/SliceTROCAL/' # INPUT w/CALIOP-TROP collocations
#pth_in_mod='/discover/nobackup/sgasso/Files/Satellite/MERRA2/'      # INPUT w/MERRA2-CAL collocations
#pth_out   ='/discover/nobackup/sgasso/Files/Satellite/Triple/'      # OUPUT w/MERRA-CAL-TROP collocations

# Name_MPL_sites , name,lat, lon, altitude(km), time1,time2
# GSFC                ,38.9930,-76.8400,0.050 , 17:05,19:05
# Santa_Cruz_Tenerife ,28.4720,-16.2470,0.052  ,13:15,15:15

############# Setup inputs
file_path  ='test_input_trop_overpass_1.txt'# File with tropomi orbit over mpl site
site_xy = [38.9930,-76.8400] # Lat/Lon MPL site
### # set up some filenames, for now make sure they are the same directory where G2GAOP is executed
config = 'm2_pm25.yaml' # this configuration file can be found in src/config
m2data = 'inst3_3d_aer_Nv' # do not change this 

### Output name for MERRA2 data along CAL track
# out_MERRA = 'm2_calipso_sampled_'+CAL_IOWA_File[12:-3]+'.nc'
save_merra='False'
out_MERRA = 'm2_test2.nc' # it may not be used depending whether to save_merra=True
### Now set output filename for triple collocation 
#output_filename = 'AlongTrack_MERRA-IOWA-CALIOP_' + caldata.cal_filename[35:-15] + '_3.1.nc'

# Get collocations dates and time of TROPO over MPL site 
data = process_file(file_path)

#####
##### now do the same for MERRA-2
#####

##### create a trajectory object
# time, lon, lat = zip(*[[data['datetimes'][0], site_xy[1], site_xy[0]],
                       # [data['datetimes'][0], site_xy[1], site_xy[0]]])
#### -----
# data_list = [[data['datetimes'][0], site_xy[1], site_xy[0]],
             # [data['datetimes'][1], site_xy[1], site_xy[0]]]
# time_list = [item[0] for item in data_list]
# lon_list  = [item[1] for item in data_list]
# lat_list  = [item[2] for item in data_list]
# time_array = np.array(time_list)
# lon_array = np.array(lon_list)
# lat_array = np.array(lat_list)

# Create time_array from the datetimes in the dictionary
time_array0 = np.array(data['datetimes'])
Nrecords   =len(time_array0)
# create time array to store window of time similar to mpl sampling
time_array =np.full((Nrecords, 3), None, dtype=object)
me_hw532   =np.zeros([Nrecords,3])
time_array[:,1]=time_array0 
# Create lon_array by repeating the first element of site_xy
lon_array = np.full(Nrecords, site_xy[1])
lat_array = np.full(Nrecords, site_xy[0])

# Define early and later arrays
delta_30_minutes = timedelta(minutes=30)
time_array[:,0] = time_array0 - delta_30_minutes
time_array[:,2] = time_array0 + delta_30_minutes

for i in [0,1,2]:
    traj = TRAJECTORY(time_array[:,i],lon_array,lat_array,m2data) # sample the MERRA-2 dataset along the trajectory, and return an xarray dataset
    #traj2 = STATION(['GSFC','GSFC'],lon_array,lat_array,m2data,time_array)
    # sample the MERRA-2 dataset along the trajectory, and return an xarray dataset
    traj_ds = traj.sample()
    print("DONE extracting MERRA data ")
    #sys.exit()
    if save_merra == 'True':
       print('\n Saving MERRA data file.... ', out_MERRA)
       traj_ds.to_netcdf(out_MERRA)
       print(" Done Savings!\n")
    #sys.exit()
    # ####### -------------------
    # # read the sampled aerosol profile data and optical tables
    optics = G2GAOP(out_MERRA,config=config)
    # Get BEXT at 532 nm and 1064nm
    # breakpoint()
    ext532  = optics.getAOPext(wavelength=532)

    # # Get MERRA-2 vertical coordinate
    ext532.pipe(addVertCoord)
    # ### Now compute MERRA weighted height
    me_zlayer =ext532.Z[0].values   # in km
    me_text532 =ext532.EXT.values # in km-1
    if len(np.where(me_text532<0)[0])>1 : sys.exit('Some EXT are <0 in MERRA. Check')
    A3=np.sum(me_zlayer * me_text532,axis=1)
    B3=np.sum(me_text532,axis=1)
    me_hw532[:,i]=np.divide(A3 , B3) # km , size=Nrecords x 1


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
