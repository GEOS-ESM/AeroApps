#!/usr/bin/env python3
"""
Nov/19/2024
This works
Test: Adaptation of calipso_sampler.ipynb notebook to AOCH applications

path to the MERRA files
/discover/nobackup/projects/gmao/merra2/data/products/MERRA2_all

"""

from pyobs.sampler import TRAJECTORY
from datetime import datetime
#from pyobs.calipso_l2 import CALIPSO_L2
import netCDF4 as nc
import numpy as np
import sys

 
from netCDF4 import Dataset,num2date
import numpy as np

class CALDATA:
    def __init__(self):
        self.filename = None
        self.lat = None
        self.lon = None
        self.time = None
        self.ext = None
        self.tback = None

def load_caldata(file_path):
    caldata = CALDATA()
    
    with Dataset(file_path, 'r') as nc:
        # Extract calipso_filename
        caldata.filename = nc.variables['calipso_filename'][0]
        caldata.lat = nc.variables['cal_lat'][:]
        caldata.lon = nc.variables['cal_lon'][:]
        caldata.ext = nc.variables['cal_ext'][:, :]
        caldata.tback = nc.variables['cal_tback'][:, :]
        # Extract and convert cal_time
        time_var = nc.variables['cal_time']
        times = time_var[:]
        time_units = time_var.units
        time_calendar = time_var.calendar
        # Convert to standard Python datetime objects
        datetime_objects = num2date(times, units=time_units, calendar=time_calendar)
        caldata.time = np.array([datetime(dt.year, dt.month, dt.day, 
                                          dt.hour, dt.minute, dt.second) 
                                 for dt in datetime_objects], dtype=object)


        #datetime_objects = num2date(times, units=time_units, calendar=time_calendar)
        #caldata.time = np.array([np.datetime64(dt) for dt in datetime_objects])
        #caldata.time = num2date(times, units=time_units, calendar=time_calendar)
    return caldata

#-----------------------------------------

# Usage
calipsoFile = 'IOWA-CALIOP_2021-08-28.nc'  
caldata = load_caldata(calipsoFile)

# Now you can access the data like this:
print("CALIPSO Filename:", caldata.filename)
print("Latitude   shape:", caldata.lat.shape)
print("Longitude  shape:", caldata.lon.shape)
print("Time       shape:", caldata.time.shape)
print("Extinction     coefficient shape:", caldata.ext.shape)
print("Backscattering coefficient shape:", caldata.tback.shape)



# set up some file names
m2data = 'inst3_3d_aer_Nv' # do not change this 
#
outFile = 'm2_calipso_sampled_'+calipsoFile[12:-3]+'.nc4'


# Read the CALIPSO file, get the lat, lon, times from the file
# caldata = CALIPSO_L2(calipsoFile,Verbose=True)


# create a trajectory object
time, lon, lat = caldata.time, caldata.lon[:], caldata.lat[:]

traj = TRAJECTORY(time,lon,lat,m2data)

# sample the MERRA-2 dataset along the trajectory, and return an xarray dataset
traj_ds = traj.sample()

# write sampled data to a netcdf file
print('Saving file.... ', outFile)
traj_ds.to_netcdf(outFile)
print('DONE!')








