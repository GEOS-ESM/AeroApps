#!/usr/bin/env python

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import yaml
from pyobs.icartt import ICARTT
import os, sys
from glob import glob
import netCDF4 as nc

#This code will generate a flight average vertical profile for a single ASIA-AQ flight. The syntax for running
#is ./plotaerosolprofiles.py VAR yyyymmdd where VAR must be a configuration listed in variablemap.yaml and the
#date must be a date in which an ASIA-AQ DC8 flight occurred as listed in plotconfig.yaml. Prior to running the
#python path needs to be set using "export PYTHONPATH=../AeroApps/install/lib/Python/" (or the full path to your
#AeroApps build).

yamlkey_var=sys.argv[1]
flightdate=sys.argv[2]

with open('variablemap.yaml') as f:
        var = yaml.safe_load(f)
with open('plotconfig.yaml') as f:
        config = yaml.safe_load(f)
        
sampledfile=glob(config['geospath'] + '/asiaaq-mrg60_dc8_*' + flightdate + '*' + var[yamlkey_var]['sampling'] + '.nc4')
ictfile=glob(config['ictpath'] + '/asiaaq-mrg60_dc8_*' + flightdate + '*ict')

# Load data
ds = xr.open_dataset(sampledfile[0], engine='netcdf4')
alt_geos = ds['H'].values
vardata_geos= ds[var[yamlkey_var]['geosvar']].values
if yamlkey_var=='OC':
        vardata_geos=vardata_geos+ds['BR'].values
ds.close
vardata_geos=vardata_geos*var[yamlkey_var]['unitconversion']

m = ICARTT(ictfile)
#print([key for key in m.__dict__.keys()])
ALT,obs_ts=m.Altitude_AGL_m_DIGANGI, m.__dict__[var[yamlkey_var]['obsvar']]

# Get the indices of minimum distances for all columns at once
indices = np.array([np.argmin(np.abs(alt_geos[n,:] - ALT[n])) for n in range(len(ALT))])

# Use advanced indexing to get the corresponding values from GEOS and filter for available obs
geos_ts = np.array([vardata_geos[n,indices[n]] for n in range(len(ALT))])
geos_ts[np.isnan(obs_ts)] = np.nan

#Find layer stats
bin_edges = np.arange(0, 4001, 250)
heightbins = np.digitize(ALT, bin_edges)
obs = np.array([np.nanmean(obs_ts[heightbins == h]) if np.any(heightbins == h) 
                         else np.nan for h in range(1, 17)])
geos = np.array([np.nanmean(geos_ts[heightbins == h]) if np.any(heightbins == h) 
                         else np.nan for h in range(1, 17)])
obs25 = np.array([np.nanpercentile(obs_ts[heightbins == h],25) if np.any(heightbins == h) 
                         else np.nan for h in range(1, 17)])
geos25 = np.array([np.nanpercentile(geos_ts[heightbins == h],25) if np.any(heightbins == h) 
                         else np.nan for h in range(1, 17)])
obs75 = np.array([np.nanpercentile(obs_ts[heightbins == h],75) if np.any(heightbins == h) 
                         else np.nan for h in range(1, 17)])
geos75 = np.array([np.nanpercentile(geos_ts[heightbins == h],75) if np.any(heightbins == h) 
                         else np.nan for h in range(1, 17)])
                         
                         
#Plot
fig = plt.figure(figsize=(10, 8))
fig.tight_layout(pad=1.0)
plt.rcParams['font.size'] = '14'
ax2 = plt.subplot(111)              
height=np.arange(0.125, 4.1, step=0.25)
plt.plot(obs,height,'black')
plt.plot(geos,height,'red')
ax2.fill_betweenx(height,obs25,obs75,facecolor='dimgray')
ax2.fill_betweenx(height,np.squeeze(geos25),np.squeeze(geos75),facecolor='red',alpha=0.3)
plt.xlabel(var[yamlkey_var]['fullname'] + ' ' + var[yamlkey_var]['units'])
plt.ylabel('Height (km)')
ax2.legend(['Obs',config['expname'],'25th-75th Percentile'],fontsize=12)
country = config['flightdate'][int(flightdate)]['country']
datestr = config['flightdate'][int(flightdate)]['datestring']
plt.title(f"{country}, {datestr}")

plt.savefig(config['expname'] + '_' + var[yamlkey_var]['geosvar'] + flightdate + '.png')
plt.show()
