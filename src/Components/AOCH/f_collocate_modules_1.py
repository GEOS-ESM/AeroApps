"""
Dec/23/2024 Moved all functions that were at the beggiing of CALTROP_ovelap_4.0.py
into this external module. So now CALTROP_ovelap_4.0.py has to read it. 


"""
import numpy as np                      # For mathematical operations (used in haversine)
import platform                         # For system and platform information (used in get_paths)
import os                              # For file path operations (used in get_paths)
import sys                             # For system exit (used in get_paths)
import re                              # For regular expressions (used in get_paths)
from netCDF4 import Dataset           # For NetCDF file operations (used in save_matched_data)
from netCDF4 import date2num         # For converting dates to numbers (used in save_matched_data)
import time                           # For timestamps in NetCDF metadata

def haversine(lat1, lon1, lat2, lon2):
        """
        Calculate the great circle distance between two points
        on the earth (specified in decimal degrees)

        lat1, lon1, lat2, lon2 are arrays
        """
        # convert decimal degrees to radians
        lon1, lat1, lon2, lat2 = map(np.radians, [lon1, lat1, lon2, lat2])

        # haversine formula
        dlon = lon2 - lon1
        dlat = lat2 - lat1
        a = np.sin(dlat/2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon/2)**2
        c = 2 * np.arcsin(np.sqrt(a))
        r = 6371 # Radius of earth in kilometers. Use 3956 for miles
        return c * r

def get_paths(yyyy, mm, dd):
    """
    Determine file paths based on operating system and computer name.
    
    Parameters:
    yyyy (str): Year
    mm (str): Month
    dd (str): Day
    
    Returns:
    tuple: (os_base_pth, tro_pthin, cal_full_path, pth_out)
    """
    
    # Initialize variables
    os_base_pth = ''
    tro_pthin = ''
    cal_full_path = ''
    pth_out = ''
    
    # Get the operating system name and computer name
    current_os = platform.system()
    computer_name = platform.node()
    print('Currently working from ', computer_name)
    
    if current_os == 'win32':  # PC office
        pthin = "D:\\Satellite\\Tropomi\\Level2\\226\\"
        
    elif current_os == 'darwin':  # Laptop
        cal_pthin = '/Volumes/ExtData1/SatData/Tropomi/Level2/'
        tro_pthin = '/Volumes/ExtData1/SatData/Tropomi/Iowa/'
        
    elif current_os == 'linux' and "calculon" in computer_name:  # Calculon
        os_base_pth = '/nobackup/CALIPSO/Level1.5/Santiago'
        cal_pthin = '/Level2/'
        tro_pthin = '/IOWA/'
        pth_out = '/nobackup/CALIPSO/Level1.5/Santiago/Output/Collocations/'
        
    elif current_os == 'Linux' and "discover" in computer_name:  # Discover
        os_base_pth = '/discover/nobackup/sgasso/Files/Satellite'
        tro_pthin = '/Iowa/'
        cal_pthin = '/css/calipso-l2-aer/data/LID_L2_05kmAPro-Standard-V4-51/'
        
        cal_folder_path = cal_pthin + yyyy + '/' + mm + '/'
        
        # Find matching files
        cal_full_path = []
        for filename in os.listdir(cal_folder_path):
            if re.search(f"{yyyy}-{mm}-{dd}..........D\.hdf$", filename):
                cal_full_path.append(os.path.join(cal_folder_path, filename))
                
        pth_out = os_base_pth + '/SliceTROCAL/'
        
    else:
        print('Current operating system not recognized.')
        print('Cannot set path to Level1 files. Terminate Run')
        sys.exit()
    
    return os_base_pth, tro_pthin, cal_full_path, pth_out

def save_matched_data(pth_out, output_file, cal_sel, tro_sel, matching_indexes, calipso_indexes, 
                     calipsol2_filename, iowa_trop_filename):
    """
    Save matched CALIOP and TROPOMI data to a NetCDF file.
    
    Parameters
    ----------
    pth_out : str
        Output directory path
    output_file : str
        Name of the output file
    cal_sel : dict
        Dictionary containing CALIOP data with keys: 'lat', 'lon', 'time', 'zlidarkm', 'ext', 'ext2', 'tback'
    tro_sel : dict
        Dictionary containing TROPOMI data with keys: 'lat', 'lon', 'height', 'aod'
    matching_indexes : array
        Array containing matching indexes between CALIOP and TROPOMI data
    calipso_indexes : array
        Array containing start and end indexes for CALIOP data
    calipsol2_filename : str
        Original CALIPSO L2 filename
    iowa_trop_filename : str
        Original TROPOMI filename
    """
    
    with Dataset(pth_out+output_file, 'w', format='NETCDF4') as nc:
        # Create dimensions
        dim_cal = nc.createDimension('dim_cal', len(cal_sel['lat']))
        dim_height = nc.createDimension('dim_height', 399)
        dim_tro = nc.createDimension('dim_tro', 2)

        # Create variables for CALIOP data
        lat_cal = nc.createVariable('cal_lat', 'f4', ('dim_cal',))
        lon_cal = nc.createVariable('cal_lon', 'f4', ('dim_cal',))
        # Create the time variable
        time_cal = nc.createVariable('cal_time', 'f8', ('dim_cal',))
        time_cal.units = 'seconds since 1970-01-01 00:00:00'
        time_cal.calendar = 'gregorian'
        time_cal.long_name = 'Time stamp in seconds since 1970-01-01'
        zlidarkm = nc.createVariable('cal_layer_km', 'f4', ('dim_height',))
        zlidarkm.long_name = "CALIPSO Layer height in km"
        ext = nc.createVariable('cal_ext', 'f4', ('dim_cal', 'dim_height'))
        ext.long_name = "Extinction Coefficient 532 nm"
        ext2 = nc.createVariable('cal_ext2', 'f4', ('dim_cal', 'dim_height'))
        ext2.long_name = "Extinction Coefficient 1064 nm"
        tback = nc.createVariable('cal_tback', 'f4', ('dim_cal', 'dim_height'))
        tback.long_name = "Backscattering Coefficient 532 nm"

        # Convert datetime objects to numeric values
        time_values = date2num(cal_sel['time'], units=time_cal.units, calendar=time_cal.calendar)

        # Assign CALIOP data
        lat_cal[:] = cal_sel['lat']
        lon_cal[:] = cal_sel['lon']
        zlidarkm[:] = cal_sel['zlidarkm']
        ext[:,:] = cal_sel['ext']
        ext2[:,:] = cal_sel['ext2']
        tback[:] = cal_sel['tback']
        time_cal[:] = time_values

        # Create variables for TROPOMI data
        lat_tro = nc.createVariable('tro_lat', 'f4', ('dim_cal',))
        lon_tro = nc.createVariable('tro_lon', 'f4', ('dim_cal',))
        height_tro = nc.createVariable('tro_height', 'f4', ('dim_cal',))
        aod_tro = nc.createVariable('tro_aod680', 'f4', ('dim_cal',))

        # Assign TROPOMI data
        lat_tro[:] = tro_sel['lat']
        lon_tro[:] = tro_sel['lon']
        height_tro[:] = tro_sel['height']
        aod_tro[:] = tro_sel['aod']

        # Create and assign string variables
        calipso_filename = nc.createVariable('calipso_filename', str)
        iowa_filename = nc.createVariable('iowa_filename', str)
        calipso_filename[0] = calipsol2_filename
        iowa_filename[0] = iowa_trop_filename

        # Create and assign array variables
        calipso_indexes_var = nc.createVariable('calipso_indexes', 'i4', ('dim_tro',))
        matching_indexes_var = nc.createVariable('matching_indexes', 'i4', ('dim_cal', 'dim_tro'))
        calipso_indexes_var[:] = np.array(calipso_indexes)
        calipso_indexes_var.description = "Start and last pixel number in selected CALIOP file that spans IOWA data"
        matching_indexes_var[:] = matching_indexes
        matching_indexes_var.description = "Array of length(CALIOP_Data) x 2 with IOWA line and column numbers that contain respective CALIOP pixels"

        # Add global attributes
        nc.title = "Combined CALIOP and TROPOMI Dataset"
        nc.description = "This file contains matched CALIOP and TROPOMI data."
        nc.history = "Created " + time.ctime(time.time())
        nc.source = "CALIOP L2 and TROPOMI L2 data products. See filenames for source files"

   # print(f"Data saved to {pth_out+output_file}")