#
# Aug/18/2025 - extract_merra_site_2.py: loops over input files.
# Jul/21/2025 - extract_merra_site_1py: based on _6, it ingests TROPO_MPL_AOCH collocations and extracts merra data. 
#_6: Jul7: modified input module so now it can read new files that include UVAI. Also cleaned some uneeded comments
#_5: July/3 = based on _4, now save to output file
#_4: Created Jul/1 = test now with actual mpl inputrs, it works
#_3: Created Jun/30 based on _2 , it computes the weigthed height 
#    This code works. Note that after making it work, I could use setup_env
# _2: Created Jun/30-
#    It tests if looping over each date, in  stn = STATION and in stn_ds = stn.sample(
#    the code runs faster. 
#    This code loops over each date of interest and it extracts MERRA data successfully
#    and fast.
# _1 Created Jun/25/2025
#    Python scripts based on mpl_sampler.ipynb
#    Also,it tests the speed of getting the data from model , particularly by changing the time
#    sampling between cases, see line 35-37

# 


import sys
####sys.path.append('/discover/nobackup/sgasso/PROJECTS/TESTJupyter/GMAOpyobs/install/lib/Python')
#sys.path.append('/gpfsm/dnb33/sgasso/REPOS/AOCH/GMAOpyobs/install/lib/Python')
#sys.path.append('/gpfsm/dnb33/sgasso/REPOS/AOCH/GMAOpyobs/src')

from pyobs.sampler import STATION
#from pyobs.mpl import MPL_L15
from pyobs.sampler import addVertCoord #for getting MERRA-2 vertical coordinate
from datetime import datetime,timedelta
import numpy as np
import time as timer
import csv,os,glob

from pyobs.aop import G2GAOP

########## Modules

# def parse_datetime(date_string):
    # try:
        # # Try with microseconds format first
        # return datetime.strptime(date_string, '%Y-%m-%d %H:%M:%S.%f')
    # except ValueError:
        # # If that fails, try without microseconds
        # return datetime.strptime(date_string, '%Y-%m-%d %H:%M:%S')

def f_parse_datetime(date_string):
    # Remove any fractional seconds that might be present
    if '.' in date_string:
        date_string = date_string.split('.')[0]
    
    # Parse using the format without fractional seconds
    return datetime.strptime(date_string, '%Y-%m-%d %H:%M:%S')
    
# def f_process_file(file_path):
    # # Read the entire file
    # with open(file_path, 'r') as file:
        # # Read all lines, skip the header
        # lines = file.readlines()[1:]

    # # Initialize lists for each column
    # datetimes = []
    # orbit_numbers = []
    # source_filenames = []
    # lines_col = []
    # columns_col = []
    # uvai_col   =[]

    # # Process each line
    # for line in lines:
        # # Split the line by comma
        # parts = line.strip().split(',')
        
        # # Append each part to its respective list
        # datetimes.append(f_parse_datetime(parts[0]))
        # orbit_numbers.append(int(parts[1]))
        # source_filenames.append(parts[2])
        # lines_col.append(int(parts[3]))
        # columns_col.append(int(parts[4]))
        # uvai_col.append(float(parts[5]))

    # return {
        # 'datetimes': datetimes,
        # 'orbit_numbers': orbit_numbers,
        # 'source_filenames': source_filenames,
        # 'lines': lines_col,
        # 'columns': columns_col,
        # 'uvai':uvai_col
    # }


def f_read_input_file(file_path):
#def read_csv_data(filename):
    """
    Read CSV file and return data as lists with datetime conversion.
    """
    try:
        datetime_list = []
        all_data = []
        
        with open(file_path, 'r') as file:
            reader = csv.reader(file)
            header = next(reader)  # Skip header
            
            for row in reader:
                # Convert date string to datetime
                dt = datetime.strptime(row[0], '%Y-%m-%d %H:%M:%S')
                datetime_list.append(dt)
                
                # Store the full row as a list
                all_data.append(row)
        
        datetime_array = np.array(datetime_list)
        
        print(f"Successfully read {len(datetime_array)} rows from {file_path}")
        return datetime_array, all_data, header
        
    except Exception as e:
        print(f"Error reading file: {e}")
        return None, None, None    
    
    
    
    
    
############# Setup inputs
#file_path  ='test_input_trop_overpass_12lines.txt'# File with tropomi orbit over mpl site
#file_path  ='GSFC_2024_TROPO_MPL_aoch.csv'# File with tropomi orbit over mpl site
#file_path  ='GSFC_2023_TROPO_MPL_aoch.csv'# File with tropomi orbit over mpl site
# file_path='/discover/nobackup/sgasso/Files/MPL/TROPO_MPL/GSFC_2018_UVAI.txt'
# file_path='/discover/nobackup/sgasso/Files/MPL/TROPO_MPL/GSFC_2024_UVAI.txt'
site_xy = [38.9930,-76.8400] # Lat/Lon MPL site

site='GSFC'

# Create a pattern to search for files containing the string
search_pattern = os.path.join('/discover/nobackup/sgasso/PROJECTS/TEST/IN/', f"*{site}*.csv")
file_list = glob.glob(search_pattern)
if not file_list : sys.exit('Stop . no input files found') 

start_time_main0 = timer.perf_counter()
# file_list=file_list[0:1]

for file_path in file_list: 
    # Extract just the filename from the full path
    filename = os.path.basename(file_path)
    
    ### Now set output filename
    #string_out='GSFC_2023'
    string_out=filename[:-19] # get input filename Site and Year information
    output_filename = 'MERRA2_TROP_MPL_site_' + string_out + '.txt'

    ### # set up some filenames, for now make sure they are the same directory where G2GAOP is executed
    config = 'm2_pm25.yaml' # this configuration file can be found in src/config
    ## Create a station object with aerosol data
    m2data = ['inst3_3d_aer_Nv','inst3_3d_asm_Nv']  # MERRA-2 collection ctrl files



    # Get collocations dates and time of TROPO over MPL site 
    # data_dict = f_process_file(file_path)
    time_array, data_list, original_header = f_read_input_file(file_path)

    # Create time_array from the datetimes in the dictionary
    #time_array = np.array(data_dict['datetimes'])


    #time_range = [time_array0.min(),time_array0.max()]
    N   =len(time_array)
    #N=2
    print(f" Read {N:6d} records \n") 
    # create time array to store window of time similar to mpl sampling
    #time_array =np.full((Nrecords, 1), None, dtype=object)
    me_hw532    =np.zeros(N)
    me_odtot532 =np.zeros(N)
    # Initialize a list of lists with N rows and 3 columns (int, float, float)
    data_output_list = [] 
    #sys.exit()

    #time_range = [time.min(),time.max()]
    #time_range = [time_array.min(),time_array.max()]
    lon, lat = site_xy[1], site_xy[0]
     
    start_time0 = timer.perf_counter()
    for i in range(N):
        time_range = [time_array[i]-timedelta(hours=1),time_array[i]+timedelta(hours=1)]
        time=np.array([time_array[i]])
        print(f"{i}  Start stn=STATION... ")
        start_time = timer.perf_counter()
        stn = STATION(['GSFC'],[lon],[lat],m2data,time_range=time_range,times=time,verbose=False)
        end_time = timer.perf_counter()
        elapsed_time = (end_time - start_time)
        print(f"           Elapsed Time= {elapsed_time:.2f} seconds")
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
        
        print("    Start stn.sample... ")
        start_time = timer.perf_counter()
        stn_ds = stn.sample(Variables=Variables).compute()
        end_time = timer.perf_counter()
        elapsed_time = (end_time - start_time)
        print(f"           Elapsed Time= {elapsed_time:.2f} seconds")
        
        ## Read the sampled aerosol profile data and optical tables
        # set up some filenames
        # this configuration file can be found in src/config
        config = 'm2_pm25.yaml'
        print("    Start G2GAOP... ")
        start_time = timer.perf_counter()
        optics = G2GAOP(stn_ds,config=config)
        end_time = timer.perf_counter()
        elapsed_time = (end_time - start_time)
        print(f"           Elapsed Time= {elapsed_time:.2f} seconds")


    # ## optional: you can write sampled data to a netcdf file
    # outFile = 'm2_mpl_sampled.nc4'
    # comp = dict(zlib=True)
    # encoding = {var: comp for var in stn_ds.data_vars}
    # stn_ds.to_netcdf(outFile,engine='netcdf4',encoding=encoding)


    # caluclate profiles of total:
    #        EXT:     aerosol extinction profile
    #        SCA:     aerosol scattering profile
    #        BSC:     aerosol backscatter profile
    #        DEPOL:   aerosol depolarization ratio
    #        ABACKTOA: total attenuated backscatter from the TOA
    #        ABACKSFC: total attenuated backscatter from the surface
    # at 532 nm
        print("    Start optics.getAOPext... ")
        start_time = timer.perf_counter()
        ext532 = optics.getAOPext(wavelength=532)
        end_time = timer.perf_counter()
        elapsed_time = (end_time - start_time)
        print(f"           Elapsed Time= {elapsed_time:.2f} seconds")


        ##### Get MERRA-2 vertical coordinate
        print(f"    Start ext532.pipe")
        start_time = timer.perf_counter()
        ext532.pipe(addVertCoord)
        end_time = timer.perf_counter()
        elapsed_time = end_time - start_time
        print(f"           Elapsed Time= {elapsed_time:.4f} seconds")
        
        ### Now compute MERRA weighted height
        print(f"    Start compute aoch")
        start_time = timer.perf_counter()
        me_zlayer =ext532.Z[0].values   # in km
        me_text532 =ext532.EXT.values # in km-1
        if len(np.where(me_text532<0)[0])>1 : sys.exit('Some EXT are <0 in MERRA. Check')
        A3=np.sum(me_zlayer[:,0] * me_text532[0,:,0],axis=0)
        B3=np.sum(me_text532[0,:,0],axis=0)
        me_hw532      =np.divide(A3 , B3) # km , size=Nrecords x 1
        me_odtot532   =np.trapz(me_text532[0,::-1,0],me_zlayer[::-1,0])
        end_time = timer.perf_counter()
        elapsed_time = end_time - start_time
        print(f"           Elapsed Time= {elapsed_time:.4f} seconds")
        
        
        
        ##### put together in output array
        #row = [data_dict['orbit_numbers'][i], me_hw532, me_odtot532]
        # data_output_list.append(row)
        row = [f"{me_hw532:.2f}",f"{me_odtot532:.2f}"]
        # Append the new values to the current row
        data_list[i].extend(row)  # or use append for individual values

        print(' ')

    # Your custom header
    #header = ['Orbit', 'AOCH532', 'TAU532']
    # Create new header with additional columns
    new_header = original_header + ['m2_aoch532', 'm2_aod532']

    # Save to CSV file
    with open(os.path.join('OUT', output_filename), 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(new_header)  # header in row 0
        writer.writerows(data_list)
    print('data saved to ', output_filename)
    elapsed_time0 = (timer.perf_counter() - start_time0)/60
    print(f"\n\n Elapsed Time= {elapsed_time0:.2f} minutes")
    print("-----------------------------")
print('Done!')

elapsed_time0 = (timer.perf_counter() - start_time_main0)/60
print(f"\n\n Total Elapsed Time= {elapsed_time0:.2f} minutes")
sys.exit()

# optional: you can write this to a netcdf file

# outFile = 'm2_mpl_ext532.nc4'
# print("Start saving to  ", outFile)
# start_time = timer.perf_counter()
# encoding = {var: comp for var in aop532.data_vars}
# aop532.to_netcdf(outFile,engine='netcdf4',encoding=encoding)
# end_time = timer.perf_counter()
# elapsed_time = (end_time - start_time)
# print(f"Saving to File: {elapsed_time:.4f} minutes")
