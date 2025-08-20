# _4: Created Jul/1 = test now with actual mpl inputrs
# _3: Created Jun/30 based on _2 , it computes the weigthed height 
#     This code works. Note that after making it work, I could use setup_env
# _2: Created Jun/30-
#     It tests if looping over each date, in  stn = STATION and in stn_ds = stn.sample(
#     the code runs faster. 
#     This code loops over each date of interest and it extracts MERRA data successfully
#     and fast.
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
from datetime import datetime
import numpy as np
import time as timer

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
    
def f_process_file(file_path):
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
        datetimes.append(f_parse_datetime(parts[0]))
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

############# Setup inputs
file_path  ='test_input_trop_overpass_12lines.txt'# File with tropomi orbit over mpl site
#file_path  ='test_input_trop_overpass_2lines.txt'# File with tropomi orbit over mpl site
#file_path='/discover/nobackup/sgasso/Files/MPL/TROPO_MPL/GSFC_2018__UVAI.txt'
site_xy = [38.9930,-76.8400] # Lat/Lon MPL site
### # set up some filenames, for now make sure they are the same directory where G2GAOP is executed
config = 'm2_pm25.yaml' # this configuration file can be found in src/config
## Create a station object with aerosol data
m2data = ['inst3_3d_aer_Nv','inst3_3d_asm_Nv']  # MERRA-2 collection ctrl files


# Get collocations dates and time of TROPO over MPL site 
data_dict = f_process_file(file_path)

# Create time_array from the datetimes in the dictionary
time_array = np.array(data_dict['datetimes'])
#time_range = [time_array0.min(),time_array0.max()]
N   =len(time_array)
# create time array to store window of time similar to mpl sampling
#time_array =np.full((Nrecords, 1), None, dtype=object)
me_hw532   =np.zeros([N,1])

#sys.exit()

## Read the MPL file
# mplAERFile = 'MPLNET_V3_L15_AER_20230605_MPL44258_GSFC.nc4'
# mplNRBFile = 'MPLNET_V3_L15_NRB_20230605_MPL44258_GSFC.nc4'
# 
# variables = ['latitude', 'longitude', 'time', 'surface_altitude', 'altitude', 'backscatter', 'c']
# mplAERdata = MPL_L15.read_file(mplAERFile,variables=variables,verb=0)
# variables = ['time','altitude','nrb']
# mplNRBdata = MPL_L15.read_file(mplNRBFile,variables=variables,verb=0)


# get the mpl obs time and location
#time       = mplAERdata.tyme
#time       = np.array([mplAERdata.tyme[0], datetime(2023, 6, 10, 1, 2, 30, 72)])
#time_array = np.array([mplAERdata.tyme[0], datetime(2023, 6, 5, 0, 2, 30, 72)])
#time_array = np.array([datetime(2023, 6, 5, 0, 0, 30, 79),datetime(2023, 6, 10, 1, 2, 30, 72)]) # two different days
#time_array = np.array([datetime(2023, 6, 5, 0, 0, 30, 79),datetime(2023, 6, 5, 2, 2, 30, 72)])  # same day, 2 hours apart
#time_array = np.array([datetime(2023, 6, 5, 0, 0, 30, 79),datetime(2023, 6, 6, 1, 2, 30, 72),datetime(2023, 6, 7, 1, 2, 30, 72)], # test 3 different dates
#time_array = np.array([datetime(2023, 6, 5, 0, 0, 30, 79),datetime(2023, 6, 6, 1, 2, 30, 72),datetime(2023, 6, 7, 1, 2, 30, 72),                       datetime(2023, 6, 8, 0, 0, 30, 79),datetime(2023, 6, 9, 1, 2, 30, 72),datetime(2023, 6,10, 1, 2, 30, 72),                       datetime(2023, 6,11, 0, 0, 30, 79),datetime(2023, 6,12, 1, 2, 30, 72),datetime(2023, 6,13, 1, 2, 30, 72)])
#time = np.array([mplAERdata.tyme[0]])

#time_range = [time.min(),time.max()]
time_range = [time_array.min(),time_array.max()]
lon, lat = site_xy[1], site_xy[0]



start_time0 = timer.perf_counter()
for i in range(N):
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
    A3=np.sum(me_zlayer * me_text532,axis=1)
    B3=np.sum(me_text532,axis=1)
    me_hw532[i]=np.divide(A3 , B3) # km , size=Nrecords x 1
    end_time = timer.perf_counter()
    elapsed_time = end_time - start_time
    print(f"           Elapsed Time= {elapsed_time:.4f} seconds")
    print(' ')




elapsed_time0 = (timer.perf_counter() - start_time0)/60
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
