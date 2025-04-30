#!/usr/bin/env python3
"""
Test code that ingests tropomi & mpl collocations, check date and location of site and pulls the MERRA data for respective
model grid cell. Then it computes AOCH  and saves in output file. 
-------------------------------------
April/30/2025 - ver 6 . It implements suggestions sent by Patricia on Apr29 to sgasso@umd.edu
April/24/2025 - ver 4 - 
-Added a loop to extract data for all years of input MPL site
-Output data is stored with TROPO and MERRA together
Apr/21/2025 - 
ver 3 works . it computes MERRA AOCH for 3 times, satellite overpass and +/- 30 min to create some sort of error bar.
But the standard deviation over the period is zero. So, I will remove this and make a version that only extracts data 
at the satellite time.

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


# def parse_datetime(date_string):
#     return datetime.strptime(date_string, '%Y-%m-%d %H:%M:%S.%f')
def parse_datetime(date_string):
    try:
        # Try with microseconds format first
        return datetime.strptime(date_string, '%Y-%m-%d %H:%M:%S.%f')
    except ValueError:
        # If that fails, try without microseconds
        return datetime.strptime(date_string, '%Y-%m-%d %H:%M:%S')
        
def process_file(file_path):
    # Read the entire file
    with open(file_path, 'r') as file:
        # Read all lines, skip the header
        lines_file = file.readlines()[1:]
    
    
    # Initialize lists for each column
    datetimes = []
    orbit_numbers = []
    source_filenames = []
    lines = []
    columns = []
    uvai = []
    ##### Process each line
    for line in lines_file:
        # Split the line by comma
        parts = line.strip().split(',')
        # Append each part to its respective list
        datetimes.append(parse_datetime(parts[0]))
        orbit_numbers.append(int(parts[1]))
        source_filenames.append(parts[2])
        lines.append(int(parts[3]))
        columns.append(int(parts[4]))
        uvai.append(float(parts[5]))

    return {
        'datetimes': datetimes,
        'orbit_numbers': orbit_numbers,
        'source_filenames': source_filenames,
        'lines': lines,
        'columns': columns,
        'uvai':uvai
    }

#-----------------------------------------

##### Set path of output file in NoNBackup so I do not fill my Discover quota
#pth_in_sat='/discover/nobackup/sgasso/Files/Satellite/SliceTROCAL/' # INPUT w/CALIOP-TROP collocations
#pth_in_mod='/discover/nobackup/sgasso/Files/Satellite/MERRA2/'      # INPUT w/MERRA2-CAL collocations
#pth_out   ='/discover/nobackup/sgasso/Files/Satellite/Triple/'      # OUPUT w/MERRA-CAL-TROP collocations

############# Setup inputs
file_path  ='test_input_trop_overpass_2.txt'# File with tropomi orbit over mpl site
#file_path='/discover/nobackup/sgasso/Files/MPL/TROPO_MPL/GSFC_2018_FullYearSites.txt'
site_xy = [38.9930,-76.8400] # Lat/Lon MPL site
### # set up some filenames, for now make sure they are the same directory where G2GAOP is executed
config = 'm2_pm25.yaml' # this configuration file can be found in src/config
m2data = 'inst3_3d_aer_Nv' # keep this file in same dir where script is exectured , do not change this 

### Output name for MERRA2 data along CAL track
# out_MERRA = 'm2_calipso_sampled_'+CAL_IOWA_File[12:-3]+'.nc'
save_merra='False'
out_MERRA = 'm2_test2.nc' # it may not be used depending whether to save_merra=True
### Now set output filename for triple collocation 
#output_filename = 'AlongTrack_MERRA-IOWA-CALIOP_' + caldata.cal_filename[35:-15] + '_3.1.nc'

# Get collocations dates and time of TROPO over MPL site 
data_dict = process_file(file_path)

#####
##### now do the same for MERRA-2
#####

# Create time_array from the datetimes in the dictionary
time_array0 = np.array(data_dict['datetimes'])
Nrecords   =len(time_array0)
# create time array to store window of time similar to mpl sampling
#time_array =np.full((Nrecords, 1), None, dtype=object)
me_hw532   =np.zeros([Nrecords,1])
#time_array[:,1]=time_array0 
# Create lon_array by repeating the first element of site_xy
lon_array = np.full(Nrecords, site_xy[1])
lat_array = np.full(Nrecords, site_xy[0])

traj = TRAJECTORY(time_array0,lon_array,lat_array,m2data) # sample the MERRA-2 dataset along the trajectory, and return an xarray dataset

### sample the MERRA-2 dataset along the trajectory, and return an xarray dataset
traj_ds = traj.sample()
print("DONE extracting MERRA data ")

### Save optical properties?
# if save_merra == 'True':
   # print('\n Saving MERRA data file.... ', out_MERRA)
   # traj_ds.to_netcdf(out_MERRA)
   # print(" Done Savings!\n")

####### -------------------
##### read the sampled aerosol profile data and optical tables
optics = G2GAOP(out_MERRA,config=config)
###Get BEXT at 532 nm and 1064nm
### breakpoint()
ext532  = optics.getAOPext(wavelength=532)
sys.exit()
##### Get MERRA-2 vertical coordinate
ext532.pipe(addVertCoord)
### Now compute MERRA weighted height
me_zlayer =ext532.Z[0].values   # in km
me_text532 =ext532.EXT.values # in km-1
if len(np.where(me_text532<0)[0])>1 : sys.exit('Some EXT are <0 in MERRA. Check')
A3=np.sum(me_zlayer * me_text532,axis=1)
B3=np.sum(me_text532,axis=1)
me_hw532=np.divide(A3 , B3) # km , size=Nrecords x 1

 
# ##### Done computing weighted heights
#### Now prepare data for saving in output file
# Create a header
header = "datetime orbit_number line column me_hw532"
# Format datetimes as strings
datetime_strings = [dt.isoformat() for dt in data_dict['datetimes']]
# Create the data array to save
data_to_save = np.column_stack((
    datetime_strings,
    data_dict['orbit_numbers'],
    data_dict['lines'],
    data_dict['columns'],
    me_hw532
))

# Save the data
filename='test_merra_collocated2.txt'
np.savetxt(filename, data_to_save, fmt='%s %5s %4s %4s %s', header=header, comments='# ')
#### use this to read data back
# Read the data back, specifying dtypes
#data = np.loadtxt(filename, dtype=[('datetime', 'U32'), ('orbit', int), ('line', int), ('column', int), ('me_hw532', float)], comments='#')

sys.exit()


# print(f"Combined height data has been added to \n {pth_out}{output_filename}")
