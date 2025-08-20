#!/usr/bin/env python3
"""
Test code that ingests tropomi & mpl collocations, check date and location of site and pulls the MERRA data for respective
model grid cell. Then it computes AOCH  and saves in output file. 
-------------------------------------
Jun26 - ver6 -- After finding a solution based on Notebook example from Patricia, this is the translation of code  test_jupyter_script_1.py into the specifics of a mpl site.


Jun/11 --- ver5 ---
Patricia created a notebook example https://github.com/GEOS-ESM/GMAOpyobs/blob/develop/src/Notebooks/sampler_examples/mpl_sampler.ipynb for sampleing over an mpl site. It uses the STATION class instead the TRAJECTORY class. 
So this ver5 implements a version of this notebook.


Apr/24/2025 - ver 4 - 
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

import sys, os
#sys.path.append('/discover/nobackup/sgasso/PROJECTS/TESTJupyter/GMAOpyobs/install/lib/Python')
sys.path.append('/gpfsm/dnb33/sgasso/REPOS/AOCH/GMAOpyobs/install/lib/Python')
sys.path.append('/gpfsm/dnb33/sgasso/REPOS/AOCH/GMAOpyobs/src')



from datetime import datetime, timedelta
import netCDF4 as nc
import numpy as np
from pyobs.sampler import STATION
from pyobs.mpl import MPL_L15
from datetime import datetime
import numpy as np
import time as timer

from pyobs.aop import G2GAOP

# import matplotlib
# import matplotlib.pyplot as plt
import shutil ### use to copy input file and then add new arrays


# def parse_datetime(date_string):
#     return datetime.strptime(date_string, '%Y-%m-%d %H:%M:%S.%f')
# def parse_datetime(date_string):
    # try:
        # # Try with microseconds format first
        # return datetime.strptime(date_string, '%Y-%m-%d %H:%M:%S.%f')
    # except ValueError:
        # # If that fails, try without microseconds
        # return datetime.strptime(date_string, '%Y-%m-%d %H:%M:%S')

def parse_datetime(date_string):
    # Remove any fractional seconds that might be present
    if '.' in date_string:
        date_string = date_string.split('.')[0]
    
    # Parse using the format without fractional seconds
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
#pth_in_sat='/discover/nobackup/sgasso/Files/Satellite/SliceTROCAL/' # INPUT w/CALIOP-TROP collocations
#pth_in_mod='/discover/nobackup/sgasso/Files/Satellite/MERRA2/'      # INPUT w/MERRA2-CAL collocations
#pth_out   ='/discover/nobackup/sgasso/Files/Satellite/Triple/'      # OUPUT w/MERRA-CAL-TROP collocations

############# Setup inputs
#file_path  ='test_input_trop_overpass_12lines.txt'# File with tropomi orbit over mpl site
file_path  ='test_input_trop_overpass_2lines.txt'# File with tropomi orbit over mpl site
#file_path='/discover/nobackup/sgasso/Files/MPL/TROPO_MPL/GSFC_2018__UVAI.txt'
site_xy = [38.9930,-76.8400] # Lat/Lon MPL site
### # set up some filenames, for now make sure they are the same directory where G2GAOP is executed
config = 'm2_pm25.yaml' # this configuration file can be found in src/config
m2data = ['inst3_3d_aer_Nv','inst3_3d_asm_Nv']  # keep this file in same dir where script is exectured , do not change this 

### Output name for MERRA2 data along CAL track
# out_MERRA = 'm2_calipso_sampled_'+CAL_IOWA_File[12:-3]+'.nc'
save_merra='False'
out_MERRA = 'm2_test2.nc' # it may not be used depending whether to save_merra=True
### Now set output filename for triple collocation 
#output_filename = 'AlongTrack_MERRA-IOWA-CALIOP_' + caldata.cal_filename[35:-15] + '_3.1.nc'

# Get collocations dates and time of TROPO over MPL site 
data_dict = process_file(file_path)
#sys.exit()
#####
##### now do the same for MERRA-2
#####

# Create time_array from the datetimes in the dictionary
time_array0 = np.array(data_dict['datetimes'])
time_range = [time_array0.min(),time_array0.max()]
Nrecords   =len(time_array0)
# create time array to store window of time similar to mpl sampling
#time_array =np.full((Nrecords, 1), None, dtype=object)
me_hw532   =np.zeros([Nrecords,1])

#time_array[:,1]=time_array0 
# Create lon_array by repeating the first element of site_xy
# lon_array = np.full(Nrecords, site_xy[1], dtype=np.float32)
# lat_array = np.full(Nrecords, site_xy[0], dtype=np.float32)


## Create a station object with aerosol data
#m2data = ['inst3_3d_aer_Nv','inst3_3d_asm_Nv']  # MERRA-2 collection ctrl files
#sys.exit()
print("Start stn.sample... ")
start_time = timer.perf_counter()
stn = STATION(['GSFC'],[site_xy[1]],[site_xy[0]],m2data,time_range=time_range,times=time_array0,verbose=True)
#stn = STATION(['GSFC'],[lon],[lat],m2data,time_range=time_range,times=time_array0,verbose=True)
end_time = timer.perf_counter()
elapsed_time = (end_time - start_time)/60
print(f"stn.sample: {elapsed_time:.4f} minutes")

#sys.exit()

## Sample the MERRA-2 dataset at the station, and return an xarray dataset
# on Jupyterhub this will take a minute - stand up, get a cup of coffee
# on the command line, this takes ~7 seconds

# Variables I want to read
du = ['DU001','DU002','DU003','DU004','DU005']
ss = ['SS001','SS002','SS003','SS004','SS005']
bc = ['BCPHILIC','BCPHOBIC']
oc = ['OCPHILIC','OCPHOBIC']
su = ['SO4']
met = ['AIRDENS','DELP','PS','RH','T']
Variables = met + du + ss + bc + oc + su

print("Start stn.sample... ")
start_time = timer.perf_counter()
stn_ds = stn.sample(Variables=Variables).compute()
end_time = timer.perf_counter()
elapsed_time = (end_time - start_time)/60
print(f"stn.sample: {elapsed_time:.4f} minutes")

print("DONE extracting MERRA data, stn.sample ")
sys.exit()
### Save optical properties?
if save_merra == 'True':
    print("Start saving into netcdf... ")
    start_time = timer.perf_counter()
    ## optional: you can write sampled data to a netcdf file
    outFile = 'm2_mpl_sampled.nc4'
    comp = dict(zlib=True)
    encoding = {var: comp for var in stn_ds.data_vars}
    stn_ds.to_netcdf(outFile,engine='netcdf4',encoding=encoding)
    end_time = timer.perf_counter()
    elapsed_time = end_time - start_time
    print(f"Done saving into netcdf, time: {elapsed_time:.4f} seconds")

    
#sys.exit()
####### -------------------
##### read the sampled aerosol profile data and optical tables
print('Now creating optics')
#optics = G2GAOP(out_MERRA,config=config)
start_time = time.perf_counter()
optics = G2GAOP(stn_ds,config=config)
end_time = time.perf_counter()
elapsed_time = end_time - start_time
print(f"G2GAOP time: {elapsed_time:.4f} seconds")

#sys.exit()
###Get BEXT at 532 nm and 1064nm
### breakpoint()
start_time = time.perf_counter()
ext532  = optics.getAOPext(wavelength=532)
end_time = time.perf_counter()
elapsed_time = end_time - start_time
print(f"optics.getAOPext time: {elapsed_time:.4f} seconds")

sys.exit()

##### Get MERRA-2 vertical coordinate
start_time = time.perf_counter()
ext532.pipe(addVertCoord)
end_time = time.perf_counter()
elapsed_time = end_time - start_time
print(f"ext532.pipe time: {elapsed_time:.4f} seconds")

### Now compute MERRA weighted height
me_zlayer =ext532.Z[0].values   # in km
me_text532 =ext532.EXT.values # in km-1
if len(np.where(me_text532<0)[0])>1 : sys.exit('Some EXT are <0 in MERRA. Check')
A3=np.sum(me_zlayer * me_text532,axis=1)
B3=np.sum(me_text532,axis=1)
me_hw532=np.divide(A3 , B3) # km , size=Nrecords x 1

sys.exit()

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
