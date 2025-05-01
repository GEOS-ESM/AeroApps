#!/usr/bin/env python3
"""
Test code that ingests tropomi & mpl collocations, check date and location of site and pulls the MERRA data for respective
model grid cell. Then it computes AOCH  and saves in output file. 
-------------------------------------
Apr/24/2025
-For user set site name , search for all matching TROPO-MPL collocations in folder , load and compute AOCH and save files
Apr/24/2025 - ver 4 - 
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
import sys, os,glob,re

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
pth_base='/discover/nobackup/sgasso/Files/MPL/'
pth_out   =pth_base + 'MERRA_MPL/'      # OUPUT 
pth_in    =pth_base + 'TROPO_MPL/'      # INPUT

### # set up some filenames, for now make sure they are the same directory where G2GAOP is executed
config = 'm2_pm25.yaml' # this configuration file can be found in src/config
m2data = 'inst3_3d_aer_Nv' # keep this file in same dir where script is exectured , do not change this 

### Output name for MERRA2 data along CAL track
# out_MERRA = 'm2_calipso_sampled_'+CAL_IOWA_File[12:-3]+'.nc'
save_merra='True'
out_MERRA = 'm2_test2.nc' # it may not be used depending whether to save_merra=True
### Now set output filename for triple collocation 
#output_filename = 'AlongTrack_MERRA-IOWA-CALIOP_' + caldata.cal_filename[35:-15] + '_3.1.nc'


############# Setup inputs
site_string='GSFC'
#file_path  ='test_input_trop_overpass_1.txt'# File with tropomi orbit over mpl site
site_xy = [38.9930,-76.8400] # Lat/Lon MPL site

# Search pattern to find all files containing the site string
search_pattern = os.path.join(pth_in, f"*{site_string}*_????_*.txt")
# Find all matching files
matching_files = glob.glob(search_pattern)
# Sort files (optional, but helps process them in a consistent order)
matching_files.sort()
print(f"Found {len(matching_files)} files for site {site_string}")

# Process each file
for input_file in matching_files:
    # Extract the base filename without path
    base_filename = os.path.basename(input_file)    # Extract year using regex pattern (looking for 4 digits after site name and underscore)
    match = re.search(f"{site_string}_(\d{{4}})_", base_filename)
    if match:
        
        # Create output filename
        # This assumes original filename has format: GSFC_2023_something.txt
        # And you want: GSFC_2023_MERRA_something.txt
        # Find what comes after the year in the original filename
        year            = match.group(1)  # Extract the year
        remaining_part  = base_filename.split(f"{site_string}_{year}_", 1)[1]
        output_filename = f"{site_string}_{year}_MERRA_{remaining_part}"
        output_path     = os.path.join(pth_out, output_filename)
        ###
        print("   Working with year " , year)
        # Get collocations dates and time of TROPO over MPL site 
        data_dict = process_file(input_file)
        
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
        print("   DONE extracting MERRA data ")
        #sys.exit()
        ### Save optical properties?
        if save_merra == 'True':
            #out_MERRA = 'm2_test2_' + site_string + str(yyyy) .nc' # it may not be used depending whether to save_merra=True
            out_MERRA = f"m2_test2_{site_string}_{year}.nc"
            print('\n Saving MERRA data file.... ', out_MERRA)
            traj_ds.to_netcdf(out_MERRA)
            print(" Done Savings!\n")
        
        ####### -------------------
        ##### read the sampled aerosol profile data and optical tables
        optics = G2GAOP(out_MERRA,config=config)
        ###Get BEXT at 532 nm and 1064nm
        ### breakpoint()
        ext532  = optics.getAOPext(wavelength=532)
        
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
        #filename='test_merra_collocated2.txt'
        np.savetxt(output_path, data_to_save, fmt='%s %5s %4s %4s %s', header=header, comments='# ')
        #### use this to read data back
        # Read the data back, specifying dtypes
        #data = np.loadtxt(filename, dtype=[('datetime', 'U32'), ('orbit', int), ('line', int), ('column', int), ('me_hw532', float)], comments='#')
    else:
        print(f"Warning: Could not extract year from filename: {base_filename}")
sys.exit()


# print(f"Combined height data has been added to \n {pth_out}{output_filename}")
