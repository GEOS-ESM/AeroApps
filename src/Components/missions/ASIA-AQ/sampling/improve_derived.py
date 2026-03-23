#!/usr/bin/env python3
"""
    Script to calcualate daily mean PM2.5 comparable to IMPROVE surface observations
    Reads in the outputs from improve_sampler.py
"""
import os, sys
import argparse
import yaml
from datetime import datetime, timedelta
import pytz
from pyobs.improve import SITE_MAP
from pyobs.xrctl import parse_ctl
from pyobs.aop import G2GAOP
import xarray as xr
import numpy as np
if __name__ == '__main__':


    parser = argparse.ArgumentParser()
    parser.add_argument("config",help='configuration yaml file')

    args = parser.parse_args()

    config = yaml.safe_load(open(args.config))

    # get IMPROVE site locations
    site_path = config['improve_site_map']
    sites = SITE_MAP(site_path)

    # get sampled data
    outdir = config['sampled_outdir'] + '/improve'
    ctl    = config['model_aer_ctl']
    sampleFile = f'{outdir}/improve.{ctl}.nc4'

    stn_ds = xr.open_dataset(sampleFile)

    # Create an optics object that links the model optics tables to the sampled aerosol profile data
    optics = G2GAOP(stn_ds,config=config['improve_config'])

    # calculate PM2.5 per species
    pm25 = {}
    for spc in optics.mieTable:
        pm25[spc] = optics.getPM(Species=[spc],pmsize=2.5,aerodynamic=True)
        # add the dry mass, we will be comparing to that later
        pmdry = pm25[spc]['PM']*(1-pm25[spc]['FWATER'])
        attrs = {'long_name':'Dry Particulate Matter Mass', 'units':'microgram m-3'}
        pmdry.attrs.update(attrs)
        pm25[spc]['PMDRY'] = pmdry

    spc = 'TOTAL'
    pm25[spc] = optics.getPM(pmsize=2.5,aerodynamic=True)
    # add the dry mass, we will be comparing to that later
    pmdry = pm25[spc]['PM']*(1-pm25[spc]['FWATER'])
    attrs = {'long_name':'Dry Particulate Matter Mass', 'units':'microgram m-3'}
    pmdry.attrs.update(attrs)
    pm25[spc]['PMDRY'] = pmdry


    # IMPROVE measures 24 hour average PM2.5 from midnight to midnight local time
    # get the daily averages at each site
    # Count valid obs per day, then mask days with fewer than 8 (3-hourly = 8/day)
    pm25_daily = {}
    for spc in pm25:
        pm25_daily[spc] = pm25[spc].sel(lev=72).resample(time='1D').mean()
        pm25_count  = pm25[spc].sel(lev=72).resample(time='1D').count()

        pm25_daily[spc] = pm25_daily[spc].where(pm25_count['PM'] >= 8)

#    sys.exit()
#    tstart = stn_ds.time[0].values.astype('datetime64[us]').item() # convert to datetime
#    tend = stn_ds.time[-1].values.astype('datetime64[us]').item()  # convert to datetime

#    tstart = tstart.replace(hour=0,minute=0)
#    tstart_iter = tstart
#    tend   = tend.replace(hour=0,minute=0)


#    # figure out how to convert local times to UTC for each site
#    localizer = []
#    for site in sites.df['SiteCode']:
#        timezone_name = sites.df[sites.df['SiteCode'] == site]['tz_name'].item()
#        if timezone_name:
#            localizer.append(pytz.timezone(timezone_name))
#        else:
#            localizer.append(None)

#    # infer time frequency of files
#    dt = xr.infer_freq(pm25['TOTAL'].time)
#    if dt == 'h':
#        ndt_per_day = 24
#    else:
#        dthours = int(dt[0])
#        ndt_per_day = 24/dthours

    
#    pm25_daily = {}
#    for spc in pm25:
#        print(spc)
#        ds = pm25[spc]
#        ds_daily = []
#        tstart_iter = tstart
#        while tstart_iter < tend:
#            print(tstart_iter)
#            ds_sites = []
#            for site, site_localizer in zip(sites.df['SiteCode'],localizer):
#                if site_localizer:
#                    time_aware_local = site_localizer.localize(tstart_iter, is_dst=False)
#                    time_utc_start = time_aware_local.astimezone(pytz.utc).replace(tzinfo=None)
#                    time_utc_end = time_utc_start + timedelta(hours=24)
#                        
#                    ds_new = ds.sel(station=site,time=slice(time_utc_start,time_utc_end),lev=72).mean(dim='time').drop_vars('lev')
#
#                    ds_new = ds_new.assign_coords(time=[tstart_iter])
#
#                    ntimes = ds.sel(station=site,time=slice(time_utc_start,time_utc_end)).sizes['time']
#
#                    if (ntimes <= ndt_per_day):
#                        # complete day not available
#                        ds_new =xr.full_like(ds_new, fill_value=np.nan)
#
#
#                    ds_sites.append(ds_new)
#
#            ds_sites = xr.concat(ds_sites,dim='station',coords="different",compat='equals')
#            ds_daily.append(ds_sites)
#
#            tstart_iter += timedelta(days=1)
#
#        pm25_daily[spc] = xr.concat(ds_daily,dim='time',coords="different",compat='equals',data_vars='all')

    # Write the speciated daily means to a netcdf file
    comp = dict(zlib=True)
    for spc in pm25:
        outFile = f'{outdir}/improve.{ctl}.pm25_{spc}.nc4'
        print('writing ',outFile)
        encoding = {var: comp for var in pm25_daily[spc].data_vars}
        pm25_daily[spc].to_netcdf(outFile,engine='netcdf4',encoding=encoding)

