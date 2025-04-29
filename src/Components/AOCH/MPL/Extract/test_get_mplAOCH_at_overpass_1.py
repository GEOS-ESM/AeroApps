#!/usr/bin/env python3
# This code ingest file with MERRA AOCH and TROPOMI orbits and then extracts the MPL data 
# for the date and time of the observations and computes MPL AOCH and then saves it. 
# Created 4/25/2025
# by @sgassoumd

from datetime import datetime,timedelta
import test_f_scripts_1 as f
import sys,os
import numpy as np


# Assuming your current script is in one folder, and the script with modules
# (let's call it 'my_module_script.py') is in a different folder:
path_to_script_directory = "/home/sgasso/REPOS/AOCH/GMAOpyobs/src/pyobs/"  # Replace with the actual path
# Add the directory to sys.path
sys.path.append(path_to_script_directory)
# Now you can import the modules by their filename (without the .py extension)
from mpl import MPL_L15
from mpl import method1
from mpl import method2



# Example usage as a script
if __name__ == "__main__":
    
    # Use a fixed filename
    filename_tro = "test_merra_collocated2.txt"
    base_path    = "/nobackup/2/MPL/L15/"
    station_str  ="GSFC"
    user_defined_minutes = 30  # window of time centered in TROPO observation overpass. 

    #### Get TROP - MPL site time collocations
    records = f.f_get_records(filename_tro)
    #### Now loop over each date and
    #### 1) verify the respective mpl direcotry exist for this date
    #### 2) open mpl file, extract data at time centered over the satellite overpass
    #### 3) store and go to next date
    N = len(records.date_time)
    for i in range(N):
        v = records.date_time[i]
        year_str = str(v.year)
        month_str = "{:02d}".format(v.month)  # Format month with leading zero if needed
        day_str = "{:02d}".format(v.day)
        #folder_path = f"/{station_str}/{year_str}/{month_str}/"
        folder_path = os.path.join(base_path, station_str, year_str, month_str)
        file_name = f"MPLNET_V3_L15_AER_{year_str}{month_str}{day_str}_MPL44258_GSFC.nc4"
        file_path = os.path.join(folder_path, file_name) # Corrected line
        if not os.path.exists(folder_path):
            print(f"Folder '{folder_path}' does not exist. Skipping to the next iteration.")
            continue  # This will jump to the next datetime object in the loop
        if not os.path.exists(file_path):
            print(f"File '{file_path}' does not exist. Skipping to the next iteration.")
            continue # This will jump to the next datetime object in the loop
            
        #print("Folder found", folder_path)
        print("File found  ", file_path)
        print(f"  Observation time {v.hour}:{v.minute}:{v.second}") 
        
        ## now for this file, ingest data and subset data within the time frame of interest. 
        # Define variables to read
        variables_to_read = ['time', 'altitude', 'extinction','extinction_err']


        # Read the data
        mpl_data = MPL_L15.read_file(file_path, variables=variables_to_read, verb=0)
        
            # ### now print some misc data directly
        # print(data.extinction.shape)
        # # Access attributes
        # print(data.time_attributes)
        # print(data.extinction_attributes.get('units'))

        #------------------------------------------------
        ###
        ### Get date from file
        # time is stored in the Julian Date format of number of days since January 1st, 4713 BC , noon UTC
        # See http://quasar.as.utexas.edu/BillInfo/JulianDatesG.html
        # and online calculator https://ssd.jpl.nasa.gov/tools/jdc/#/cd
        # where integer par is Ndays since this date and decimal time is hh:mm:ss since noon UTC
        # so for example : 2440588.000000 is A.D. 1970 January 1    12:00:00.0
        # this number can be substracted from data.time and use it as new reference point
        # altenrnativel the python module astropy can deal with it. NOTE the python module
        # datetime does not have a way to deal with this reference point .
        # I used a different approach based on the date in the filename.

        x     =mpl_data.time
        base_date = datetime(v.year, v.month, v.day,0,0,0) ## Your base date
        xtime = []
        for jd in x:
            # Get just the fractional part since we know the date
            day_fraction = jd % 1
            if day_fraction >= 0.5:
                day_fraction = day_fraction - 0.5
            else:
                day_fraction = day_fraction+0.5
            # Convert fractional day to minutes
            #seconds = int(day_fraction * 86400)
            minutes = int(day_fraction * 1440)
            # Create datetime by adding the time to base date
            dt = base_date + timedelta(minutes=minutes)
            xtime.append(dt)
        ### Because MPL mpl_data is reported every 60 min,
        ### one could create a time array of size 60*24 = 1440 and it should be the same.
        # Convert the list to a numpy array
        xtime = np.array(xtime)
        
        ### Assign inputs to module
        #site_name = station_str
        z_in      = mpl_data.altitude[0,:]*1000 # convert to meters

        ### Now select variable of interest and transpose so it is compatible
        ### with plotting routine
        varin       =1e-3*mpl_data.extinction[0,:].T # 1440 x 400 to 400 x 1440
        varin_err   =1e-3*mpl_data.extinction_err[0,:].T # 1440 x 400 to 400 x 1440


        ##### Set time frame of satellite overpass
        # Create a timedelta object representing X minutes
        delta = timedelta(minutes=user_defined_minutes)
        # Calculate v + X minutes
        t_end = v + delta
        # Calculate v - X minutes
        t_start = v - delta

        # Find the indexes where the values span the time frame
        indexes = np.where((xtime >= t_start) & (xtime <= t_end))[0]

        # Extract the respective values using the indexes
        var_sub = varin[:,indexes]
        var_sub_err=varin_err[:,indexes]

        ### now proceed to compute AOCH
        # print("Indexes:", indexes)


        # Method 1 : average over time each layer and the integrate over column
        # Method 2 : Integrate each column and then take average AOCHs
        # Apply both methods
        # result1, error1 = method1(varin, varin_err, z_in)
        # result2, error2 = method2(varin, varin_err, z_in)

        # Apply both methods
        result1, error1, total_obs1, valid_obs1, used_in_integral1 = method1(var_sub, var_sub_err, z_in)
        result2, error2, total_obs2, valid_obs2, used_in_integral2 = method2(var_sub, var_sub_err, z_in)

        # Print results for Method 1
        print("  Method 1 : Average Bext in each layer and then compute AOCH")
        if np.isnan(result1):
            print("  AOCH (km): No valid MPL data for integration")
        else:
            print(f"  AOCH (km): {result1/1000:.2f} ± {error1/1000:.2f}")
            print(f"    Total Profiles: {total_obs1}")
            print(f"    Stats of N layers with EXT obs: min={valid_obs1.min()}, max={valid_obs1.max()}, mean={valid_obs1.mean():.2f}")
            print(f"    N layers(out of 400) used in final mean: {used_in_integral1}")

        print("\n  Method 2: Compute AOCH in each profiles and then take mean:")
        if np.isnan(result2):
            print("  AOCH (km): No valid MPL data for integration")
        else:
            print(f"  AOCH (km): {result2/1000:.2f} ± {error2/1000:.2f}")
            print(f"    Total Profile: {total_obs2}")
            print(f"    Stats of N layers w/EXT obs per profile: min={valid_obs2.min()}, max={valid_obs2.max()}, mean={valid_obs2.mean():.2f}")
            print(f"    N Profiles used in final value : {used_in_integral2}")
        
        print(" ")



