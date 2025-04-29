#!/usr/bin/env python3
# This code ingest file with MERRA AOCH and TROPOMI orbits and then extracts the MPL data 
# for the date and time of the observations and computes MPL AOCH and then saves it. 
# Created 4/25/2025
# 4/29/2025 - _2.py : Maine difference with _1.py is that now it saves data in output file (only for cases with MPL AOCH computed) 
#   and the main iteration is improved after chatGSFC
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
    # Create lists to store all the results
    all_results = []

    # Rerun the loop to store data
    for i in range(N):
        v = records.date_time[i]
        year_str = str(v.year)
        month_str = "{:02d}".format(v.month)
        day_str = "{:02d}".format(v.day)
        folder_path = os.path.join(base_path, station_str, year_str, month_str)
        file_name = f"MPLNET_V3_L15_AER_{year_str}{month_str}{day_str}_MPL44258_GSFC.nc4"
        file_path = os.path.join(folder_path, file_name)
        
        # Initialize result variables with NaN
        result1 = np.nan
        error1 = np.nan
        total_obs1 = np.nan
        valid_obs1 = np.array([np.nan])  # Initialize as array with NaN
        valid_obs1_min = np.nan
        valid_obs1_max = np.nan
        valid_obs1_mean = np.nan
        used_in_integral1 = np.nan
        
        result2 = np.nan
        error2 = np.nan
        total_obs2 = np.nan
        valid_obs2 = np.array([np.nan])  # Initialize as array with NaN
        valid_obs2_min = np.nan
        valid_obs2_max = np.nan
        valid_obs2_mean = np.nan
        used_in_integral2 = np.nan
        
        has_valid_data = False
           # Check if file and folder exist
        if os.path.exists(folder_path) and os.path.exists(file_path):
            # Read the data
            variables_to_read = ['time', 'altitude', 'extinction', 'extinction_err']
            mpl_data = MPL_L15.read_file(file_path, variables=variables_to_read, verb=0)
            
            
            # ### now print some misc data directly
            # print(data.extinction.shape)
            #### # Access attributes
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

            x = mpl_data.time
            base_date = datetime(v.year, v.month, v.day, 0, 0, 0)
            xtime = []
            for jd in x:
                day_fraction = jd % 1
                if day_fraction >= 0.5:
                    day_fraction = day_fraction - 0.5
                else:
                    day_fraction = day_fraction + 0.5
                minutes = int(day_fraction * 1440)
                dt = base_date + timedelta(minutes=minutes)
                xtime.append(dt)
            xtime = np.array(xtime)
            
            z_in = mpl_data.altitude[0,:] * 1000
            varin = 1e-3 * mpl_data.extinction[0,:].T
            varin_err = 1e-3 * mpl_data.extinction_err[0,:].T
            
            delta = timedelta(minutes=user_defined_minutes)
            t_end = v + delta
            t_start = v - delta
            
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

            ### Check if result1 or result2 contains non-NaN values
            valid_data_1 = not np.isnan(result1)
            valid_data_2 = not np.isnan(result2)
            ### Set flag 
            has_valid_data = valid_data_1 or valid_data_2
            
            if valid_data_1:
                valid_obs1_min = valid_obs1.min()
                valid_obs1_max = valid_obs1.max()
                valid_obs1_mean = valid_obs1.mean()
            
            if valid_data_2:
                valid_obs2_min = valid_obs2.min()
                valid_obs2_max = valid_obs2.max()
                valid_obs2_mean = valid_obs2.mean()

            # Print results for Method 1
            print("\n  Method 1 : Average Bext in each layer and then compute AOCH")
            if np.isnan(result1):
                print("    AOCH (km): No valid data for integration")
            else:
                print(f"  AOCH (km): {result1/1000:.2f} ± {error1/1000:.2f}")
                print(f"    Total Profiles: {total_obs1}")
                print(f"    Stats of N layers with EXT obs: min={valid_obs1.min()}, max={valid_obs1.max()}, mean={valid_obs1.mean():.2f}")
                print(f"    N layers(out of 400) used in final mean: {used_in_integral1}")

            print("  Method 2: Compute AOCH in each profiles and then take mean:")
            if np.isnan(result2):
                print("    AOCH (km): No valid data for integration")
            else:
                print(f"  AOCH (km): {result2/1000:.2f} ± {error2/1000:.2f}")
                print(f"    Total Profile: {total_obs2}")
                print(f"    Stats of N layers w/EXT obs per profile: min={valid_obs2.min()}, max={valid_obs2.max()}, mean={valid_obs2.mean():.2f}")
                print(f"    N Profiles used in final value : {used_in_integral2}")
        # Only create and append the row if we have valid data
        if has_valid_data:
            # Create a row with all the data you want to save
            # Start with the data from the original record
            row = {}
            for key in records.__dict__.keys():
                if hasattr(records, key) and i < len(getattr(records, key)):
                    row[key] = getattr(records, key)[i]
            
            # Add computed results
            row['method1_aoch_m'] = result1
            row['method1_error_m'] = error1
            row['method1_total_obs'] = total_obs1
            row['method1_valid_obs_min'] = valid_obs1_min
            row['method1_valid_obs_max'] = valid_obs1_max
            row['method1_valid_obs_mean'] = valid_obs1_mean
            row['method1_used_in_integral'] = used_in_integral1
            
            row['method2_aoch_m'] = result2
            row['method2_error_m'] = error2
            row['method2_total_obs'] = total_obs2
            row['method2_valid_obs_min'] = valid_obs2_min
            row['method2_valid_obs_max'] = valid_obs2_max
            row['method2_valid_obs_mean'] = valid_obs2_mean
            row['method2_used_in_integral'] = used_in_integral2
            
            all_results.append(row)
            
            
            
            print(" ")
            
    # Only proceed if we have results to save
    if all_results:
        # Determine the output filename by adding "_with_mpl" before the extension
        output_filename = os.path.splitext(filename_tro)[0] + "_with_mpl" + os.path.splitext(filename_tro)[1]

        # Create header based on the keys of the first result
        header = list(all_results[0].keys())

        # Write the data to the output file
        with open(output_filename, 'w') as f:
            # Write the header
            f.write('\t'.join(header) + '\n')
            
            # Write each row of data
            for row in all_results:
                line = []
                for key in header:
                    value = row[key]
                    
                    # Format the value based on its type
                    if isinstance(value, datetime):
                        formatted_value = value.strftime("%Y-%m-%d %H:%M:%S")
                    elif isinstance(value, (float, np.float32, np.float64)):
                        if np.isnan(value):
                            formatted_value = "NaN"
                        else:
                            formatted_value = f"{value:.6f}"
                    else:
                        formatted_value = str(value)
                        
                    line.append(formatted_value)
                
                f.write('\t'.join(line) + '\n')

        print(f"\nData saved to {output_filename} with {len(all_results)} valid records")
    else:
        print("\nNo valid MPL data found, output file was not created")



