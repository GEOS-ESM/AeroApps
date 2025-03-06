
### test_f_find_orbits_2.py: Created 2/26/2025 - 
###       Using _1 as template, this code implements a simple way 
###       of estimatating what orbits are likely to contain 
#--------------------------------------
### created Feb/12/2025 : these are the modules used by test_create_list_orbits_with_sites_2.py


import os,sys,platform
import glob
#import h5py as h5
import numpy as np
from pathlib import Path
from datetime import datetime, timedelta


def find_hdf5_files(folder_path):
    ##### Find all HDF5 files in the specified folder
    hdf5_files = []
    for file in Path(folder_path).glob('*.hdf5'):  # Also add '*.h5' if needed
        hdf5_files.append(str(file))
    return hdf5_files

def get_path(now_os, now_computer):
    #### Module to setup path to where the datafiles are. Needed if working in multiple computers
    # Your existing get_path function remains the same
    if now_os=='win32':
        base_path    = 'C:/Users/sgasso/Downloads/'
        pth_fig_out='C:/Users/sgasso/OneDrive - NASA/ToShare/2025/GEOS/Pyfigures/'
    elif now_os == 'darwin':
        base_path= '/Volumes/ExtData1/SatData/Tropomi/Level2/2023/359/'
        pth_fig_out='/Users/sgasso/Library/CloudStorage/OneDrive-NASA/ToShare/2025/AOCH/PyFigures/'
    elif now_os == 'linux' and "calculon" in now_computer:
        base_path = '/nobackup/TROPOMAER/'
        pth_fig_out = ''
    elif now_os == 'linux' and "discover" in now_computer:
        base_path = '/nobackup/CALIPSO/Level1.5/Santiago/'
    else:
        print('Current operating system no recognized.')
        print('Cannot set path to MPL  files. Terminate Run')
        sys.exit()
    return base_path,pth_fig_out


#### Read list of MPL sites
def read_mpl_site_list(filename):
#Expected format: comma separated
#line 1 :Name_MPL_sites , name,lat, lon, altiude(km),time_start, time_end
#line 2 : GSFC           ,38.9930,-76.8400,0.050,    hh:mm , hh:mm
#.............
# Last two columns are ignored 
    sites = []
    with open(filename, 'r') as file:
        header = file.readline()
        for line in file:
            data = [item.strip() for item in line.split(',')]
            site = [
                data[0],
                float(data[1]),
                float(data[2]),
                float(data[3]),
            ]
            sites.append(site)
    return sites


def get_orbit_site_matches2(list_mpl_sites,yyyy, mm=None, dd=None, julian=None, verb=True):
#--------------------------------------------------------------------------------------------
#     Module to match orbits with MPL sites based on temporal overlap.
#     
#     Parameters:
#     -----------
#     yyyy : int , Year
#     mm : int, optional ,  Month (1-12)
#     dd : int, optional ,  Day (1-31)
#     julian : int, optional , Julian day (alternative to mm/dd)
#     verb : bool, optional ,  Verbose output flag (default: True)
#         
#     Returns:
#     --------
#     list
#         List containing files with files in selected directory and the path to the folder
#--------------------------------------------------------------------------------------------
    
    if not yyyy:
        raise ValueError('Input year must be provided')
    if not mm and not julian:
        raise ValueError('Either Month and day or Julian day must be provided')

    # Set correct paths according to the current computer and OS
    current_working_directory = os.getcwd()
    current_os           = platform.system()
    computer_name        = platform.node()
    current_pth, pth_fig = get_path(current_os.lower(), computer_name)
    
    # Convert mm/dd to julian if needed
    if mm:
        date_obj = datetime(yyyy, mm, dd)
        julian = (date_obj - datetime(yyyy, 1, 1)).days + 1
    julian_str = f"{julian:03d}"

    if verb:print(f'    Year is {yyyy}, julian {julian_str}')

    # Get list of files for selected date
    path_2_folder_with_orbits = current_pth+str(yyyy)+'/'+julian_str +'/'
    file_list = glob.glob(path_2_folder_with_orbits+'*.nc')
    
    # Process orbit information
    orbits_list = []
    for full_pathname in file_list:
        orbits_list.append(full_pathname[-80:])

    return orbits_list,path_2_folder_with_orbits

