# Created Jun/25/2025
# Python scripts based on mpl_sampler.ipynb
# Also,it tests the speed of getting the data from model , particularly by changing the time
# sampling between cases, see line 35-37

# 


import sys
#sys.path.append('/discover/nobackup/sgasso/PROJECTS/TESTJupyter/GMAOpyobs/install/lib/Python')
sys.path.append('/gpfsm/dnb33/sgasso/REPOS/AOCH/GMAOpyobs/install/lib/Python')
sys.path.append('/gpfsm/dnb33/sgasso/REPOS/AOCH/GMAOpyobs/src')



from pyobs.sampler import STATION
from pyobs.mpl import MPL_L15
from datetime import datetime
import numpy as np
import time as timer

from pyobs.aop import G2GAOP

## Read the MPL file
mplAERFile = 'MPLNET_V3_L15_AER_20230605_MPL44258_GSFC.nc4'
mplNRBFile = 'MPLNET_V3_L15_NRB_20230605_MPL44258_GSFC.nc4'

variables = ['latitude', 'longitude', 'time', 'surface_altitude', 'altitude', 'backscatter', 'c']
mplAERdata = MPL_L15.read_file(mplAERFile,variables=variables,verb=0)
variables = ['time','altitude','nrb']
mplNRBdata = MPL_L15.read_file(mplNRBFile,variables=variables,verb=0)
print(mplNRBdata)

# get the mpl obs time and location
#time = mplAERdata.tyme
time = np.array([mplAERdata.tyme[0], datetime(2023, 6, 10, 1, 2, 30, 72)])# 5 days apart
#time = np.array([mplAERdata.tyme[0], datetime(2023, 6, 5, 0, 2, 30, 72)]) # 2 hours apart 
#time = np.array([mplAERdata.tyme[0]])

time_range = [time.min(),time.max()]
lon, lat = mplAERdata.longitude[0], mplAERdata.latitude[0]

## Create a station object with aerosol data
m2data = ['inst3_3d_aer_Nv','inst3_3d_asm_Nv']  # MERRA-2 collection ctrl files

print("Start stn=STATION... ")
start_time = timer.perf_counter()
stn = STATION(['GSFC'],[lon],[lat],m2data,time_range=time_range,times=time,verbose=True)
end_time = timer.perf_counter()
elapsed_time = (end_time - start_time)/60
print(f"stn = STATION([: {elapsed_time:.4f} minutes")

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

sys.exit()
## optional: you can write sampled data to a netcdf file
outFile = 'm2_mpl_sampled.nc4'
comp = dict(zlib=True)
encoding = {var: comp for var in stn_ds.data_vars}
stn_ds.to_netcdf(outFile,engine='netcdf4',encoding=encoding)


## Read the sampled aerosol profile data and optical tables
# set up some filenames
# this configuration file can be found in src/config
config = 'm2_pm25.yaml'
print("Start G2GAOP... ")
start_time = timer.perf_counter()
optics = G2GAOP(stn_ds,config=config)
end_time = timer.perf_counter()
elapsed_time = (end_time - start_time)
print(f"G2GAOP(stn_ds: {elapsed_time:.4f} minutes")

# caluclate profiles of total:
#        EXT:     aerosol extinction profile
#        SCA:     aerosol scattering profile
#        BSC:     aerosol backscatter profile
#        DEPOL:   aerosol depolarization ratio
#        ABACKTOA: total attenuated backscatter from the TOA
#        ABACKSFC: total attenuated backscatter from the surface
# at 532 nm
print("Start optics.getAOPext... ")
start_time = timer.perf_counter()
aop532 = optics.getAOPext(wavelength=532)
end_time = timer.perf_counter()
elapsed_time = (end_time - start_time)
print(f"optics.getAOPext: {elapsed_time:.4f} minutes")

# optional: you can write this to a netcdf file

outFile = 'm2_mpl_ext532.nc4'
print("Start saving to  ", outFile)
start_time = timer.perf_counter()
encoding = {var: comp for var in aop532.data_vars}
aop532.to_netcdf(outFile,engine='netcdf4',encoding=encoding)
end_time = timer.perf_counter()
elapsed_time = (end_time - start_time)
print(f"Saving to File: {elapsed_time:.4f} minutes")
