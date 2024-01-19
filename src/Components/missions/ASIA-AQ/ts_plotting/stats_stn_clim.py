import pandas as pd
import xarray as xr
from pyobs.sampler import STATION
from glob import glob
from datetime import datetime


if __name__ == "__main__":

    stations = 'seoul'
#    stations = 'manila'

    dataset = sorted(glob('{}_2*.nc'.format(stations)))
    ds = xr.open_mfdataset(dataset,parallel=True)
    ds['PM25'] = ds['BCSMASS']+ds['OCSMASS']+ds['DUSMASS25']+ds['SO4SMASS']+ds['SSSMASS25']
    ds['CCEXTTAU'] = ds['OCEXTTAU']+ds['BCEXTTAU']


    # get hourly climatology statistics
    grouper = xr.DataArray(
        pd.MultiIndex.from_arrays([ds.time.dt.month.values, ds.time.dt.day.values,ds.time.dt.hour.values],
        names=['month', 'day','hour']),
        dims=['time'],coords=[ds.time])

    clim_mean = {}
    clim_med  = {}
    clim_25   = {}
    clim_75   = {}
    clim_min  = {}
    clim_max  = {}
    for var in list(ds.keys()):
        clim_mean[var] = ds[var].groupby(grouper).mean()
        clim_min[var]  = ds[var].groupby(grouper).min()
        clim_max[var]  = ds[var].groupby(grouper).max()
        clim_25[var]   = ds[var].chunk(dict(time=-1)).groupby(grouper).quantile(0.25)
        clim_med[var]  = ds[var].chunk(dict(time=-1)).groupby(grouper).quantile(0.50)
        clim_75[var]   = ds[var].chunk(dict(time=-1)).groupby(grouper).quantile(0.75)

    clim_mean = xr.Dataset(clim_mean).rename({'group': 'time'})
    clim_25   = xr.Dataset(clim_25).rename({'group': 'time'})
    clim_75   = xr.Dataset(clim_75).rename({'group': 'time'})
    clim_med  = xr.Dataset(clim_med).rename({'group': 'time'})
    clim_min  = xr.Dataset(clim_min).rename({'group': 'time'})
    clim_max  = xr.Dataset(clim_max).rename({'group': 'time'})
    
    # create generic time coordinate
    time = [datetime(2024,m,d,h) for m,d,h in zip(clim_mean.time.time_level_0.values,clim_mean.time.time_level_1.values,clim_mean.time.time_level_2.values)]
    clim_mean = clim_mean.assign_coords({'time':time}) 
    clim_25 = clim_25.assign_coords({'time':time})
    clim_75 = clim_75.assign_coords({'time':time})
    clim_med = clim_med.assign_coords({'time':time})
    clim_min = clim_min.assign_coords({'time':time})
    clim_max = clim_max.assign_coords({'time':time})

    # write to netcdf file
    clim_mean.to_netcdf(path='{}_clim_mean.nc'.format(stations))
    clim_25.to_netcdf(path='{}_clim_25.nc'.format(stations))
    clim_75.to_netcdf(path='{}_clim_75.nc'.format(stations))
    clim_med.to_netcdf(path='{}_clim_med.nc'.format(stations))
    clim_min.to_netcdf(path='{}_clim_min.nc'.format(stations))
    clim_max.to_netcdf(path='{}_clim_max.nc'.format(stations))
