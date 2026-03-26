#!/usr/bin/env python3
"""
    Script to create plots comparing a GEOS simulations to each other sampled at IMPROVE sites 
    Reads in the outputs from sample/improve_derived.py
"""
import os, sys
import argparse
import yaml
from pyobs.improve import SITE_MAP, IMPROVE
import xarray as xr
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats
from map_plot_mae import plot_station_map_mae
from map_plot_mb import plot_station_map_mb

def clean_node(ds):
    if ds is None:
        return ds
    # Add .copy() to get a mutable version
    ds = ds.copy()
    
    for var in list(ds.variables):
        if ds[var].dtype == object:
            print(f'  Converting {var} to string')
            try:
                # Try converting to string
                ds[var] = ds[var].astype(str)
            except Exception as e:
                # If conversion fails, drop the variable
                print(f'  Could not convert {var}, dropping: {e}')
                ds = ds.drop_vars(var)
    
    return ds

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("config",help='configuration yaml file')

    args = parser.parse_args()

    config = yaml.safe_load(open(args.config))

    # get model sampled data 
    outdir = config['sampled_outdir']
    baseline_ctl    = config['model_baseline']
    oxh_ctl         = config['model_oxh']
    oxm_ctl         = config['model_oxm']
    spcs = ['BC','BR','DU','NH4','NI','OC','SS','SU','TOTAL']

    # get output directory
    plot_outdir = config['plot_outdir'] + '/' + baseline_ctl
    os.makedirs(plot_outdir,exist_ok=True)

    if config.get('do_preprocess', True):
        baseline_ds = {}
        oxh_ds = {}
        oxm_ds = {}
        for spc in spcs:
            baseline_ds[spc] = xr.open_dataset(f'{outdir}/improve.{baseline_ctl}.pm25_{spc}.nc4')
            oxh_ds[spc] = xr.open_dataset(f'{outdir}/improve.{oxh_ctl}.pm25_{spc}.nc4')
            oxm_ds[spc] = xr.open_dataset(f'{outdir}/improve.{oxm_ctl}.pm25_{spc}.nc4')

        # get time span of data
        syear = pd.Timestamp(baseline_ds['TOTAL'].time[0].values).year
        eyear = pd.Timestamp(baseline_ds['TOTAL'].time[-1].values).year

        yearlist = list(range(syear,eyear+1))

        # get IMPROVE site locations
        site_path = config['improve_site_map']
        sites = SITE_MAP(site_path)

        # get IMPROVE data
        imp_path = config['improve_data_dir']
        implist = []
        for year in yearlist:
            implist += [f'{imp_path}/IMPAER_{year}.txt']
        imp = IMPROVE(implist,site_path=site_path)

        # convert IMPROVE dataframe to xarray dataset for easier comparisons
        df = imp.df.rename(columns={'SiteCode':'station','Start_Time':'time','Longitude':'lon','Latitude':'lat'})
        df = df.set_index(['time','station','ParamCode'])
        print(f"Original length: {len(df.index)}")
        print(f"Unique index length: {len(df.index.unique())}")
        # there  can be duplicate entries because there are multiple instruments at one site
        # let's just keep the first measurement in the database
        df_unique = df[~df.index.duplicated(keep='first')]
        # Now convert the cleaned DataFrame to xarray
        imp_ds = df_unique.to_xarray()
        # set all -999 to nan
        imp_ds = imp_ds.where(imp_ds != -999.)
        # retain all pm25 values to where Status == 'V0', set all the rest to nan
        imp_ds["pm25"] = imp_ds["pm25"].where(imp_ds["Status"] == "V0")


        # subset geos and improve for the dates that are available in both
        for spc in spcs:
            # Find the common times between both datasets
            common_times = baseline_ds[spc].time.values[np.isin(baseline_ds[spc].time.values, imp_ds.time.values)]


            # align coordinates
            baseline_sync = baseline_ds[spc].sel(time=common_times,station=imp_ds.station.values)
            baseline_ds[spc] = baseline_sync

            oxh_sync = oxh_ds[spc].sel(time=common_times,station=imp_ds.station.values)
            oxh_ds[spc] = oxh_sync

            oxm_sync = oxm_ds[spc].sel(time=common_times,station=imp_ds.station.values)
            oxm_ds[spc] = oxm_sync


        baseline_ds['improve'] = imp_ds
        oxh_ds['improve'] = imp_ds
        oxm_ds['improve'] = imp_ds

        # create datatrees for each experiment
        dt_baseline = xr.DataTree.from_dict(baseline_ds)
        dt_oxh = xr.DataTree.from_dict(oxh_ds)
        dt_oxm = xr.DataTree.from_dict(oxm_ds)

        # save to zarr
        dt_clean = dt_baseline.map_over_datasets(clean_node)
        dt_clean.to_zarr(f'{outdir}/dt_baseline.zarr', mode='w')
        dt_clean = dt_oxh.map_over_datasets(clean_node)
        dt_clean.to_zarr(f'{outdir}/dt_oxh.zarr', mode='w')
        dt_clean = dt_oxm.map_over_datasets(clean_node)
        dt_clean.to_zarr(f'{outdir}/dt_oxm.zarr', mode='w')

    else:
        # get IMPROVE site locations
        site_path = config['improve_site_map']
        sites = SITE_MAP(site_path)

        dt_baseline = xr.open_datatree(f'{outdir}/dt_baseline.zarr', engine='zarr') 
        dt_oxh = xr.open_datatree(f'{outdir}/dt_oxh.zarr', engine='zarr')
        dt_oxm = xr.open_datatree(f'{outdir}/dt_oxm.zarr', engine='zarr')


    ## TOTAL PM2.5
    spc_mapping = {
        'RCFM': 'TOTAL',
        'SeaSaltf': 'SS',
        'SOILf': 'DU',
        'ECf': 'BC',
        'OMCf': 'OC',  # will need special handling for OC+BR
        'ammSO4f': 'SU',
        'ammNO3f': 'NI',
    }

#    spc_mapping = {'ammNO3f': 'NI'}

    for spc in spc_mapping:
        # Get baseline data
        if spc == 'OMCf':
            baseline_spc = dt_baseline['OC'] + dt_baseline['BR']
            oxh_spc = dt_oxh['OC'] + dt_oxh['BR']
            oxm_spc = dt_oxm['OC'] + dt_oxm['BR']
        else:
            baseline_spc = dt_baseline[spc_mapping[spc]]
            oxh_spc = dt_oxh[spc_mapping[spc]]
            oxm_spc = dt_oxm[spc_mapping[spc]]
        
        # Get observation data
        imp_data = dt_baseline['improve'].pm25.sel(ParamCode=spc, drop=True)
        
        # Create dataframes for baseline vs obs (top subplot)
        ds_baseline = baseline_spc.PMDRY.where(imp_data.notnull())
        monthly_mean = ds_baseline.groupby('time.month').mean(skipna=True).rename({'month': 'time'})
        df_baseline = ds_baseline.to_dataframe(name="PM25").reset_index()
        df_baseline["Source"] = "Baseline"
        
        monthly_mean = imp_data.groupby('time.month').mean(skipna=True).rename({'month': 'time'})
        df_imp = imp_data.to_dataframe(name="PM25").reset_index()
        df_imp["Source"] = 'IMPROVE'
        
        # Create dataframe for top subplot (baseline vs obs)
        plot_df_top = pd.concat([df_baseline, df_imp], ignore_index=True)
        
        # Create dataframes for bottom subplot (experiment differences)
        ds_oxh = oxh_spc.PMDRY.where(imp_data.notnull())
        ds_oxm = oxm_spc.PMDRY.where(imp_data.notnull())
        
        
        # Calculate differences (experiment - improve)
        ds_baseline_diff = ds_baseline - imp_data
        monthly_mean = ds_baseline_diff.groupby('time.month').mean(skipna=True).rename({'month': 'time'})
        df_baseline_diff = ds_baseline_diff.to_dataframe(name="PM25").reset_index()
        df_baseline_diff['Source'] = 'Baseline - IMPROVE'
        # Merge lat/lon from sites.df into df_baseline_diff
        # matching on SiteCode (sites.df) = station (df_baseline_diff)
        df_baseline_diff = df_baseline_diff.merge(
            sites.df[['SiteCode', 'Latitude', 'Longitude']],
            left_on='station',
            right_on='SiteCode',
            how='left'
        ).drop(columns='SiteCode')  # drop the redundant SiteCode column


        ds_oxh_diff = ds_oxh - imp_data
        monthly_mean = ds_oxh_diff.groupby('time.month').mean(skipna=True).rename({'month': 'time'})
        df_oxh_diff= ds_oxh_diff.to_dataframe(name="PM25").reset_index()
        df_oxh_diff['Source'] = 'OxH - IMPROVE'
        # Merge lat/lon from sites.df into df_baseline_diff
        # matching on SiteCode (sites.df) = station (df_baseline_diff)
        df_oxh_diff = df_oxh_diff.merge(
            sites.df[['SiteCode', 'Latitude', 'Longitude']],
            left_on='station',
            right_on='SiteCode',
            how='left'
        ).drop(columns='SiteCode')  # drop the redundant SiteCode column
        
        ds_oxm_diff = ds_oxm - imp_data
        monthly_mean = ds_oxm_diff.groupby('time.month').mean(skipna=True).rename({'month': 'time'})
        df_oxm_diff = ds_oxm_diff.to_dataframe(name="PM25").reset_index()
        df_oxm_diff['Source'] = 'OxM - IMPROVE'
        # Merge lat/lon from sites.df into df_baseline_diff
        # matching on SiteCode (sites.df) = station (df_baseline_diff)
        df_oxm_diff = df_oxm_diff.merge(
            sites.df[['SiteCode', 'Latitude', 'Longitude']],
            left_on='station',
            right_on='SiteCode',
            how='left'
        ).drop(columns='SiteCode')  # drop the redundant SiteCode column

        
        # Create dataframe for bottom subplot
        df_list = [df_baseline_diff, df_oxm_diff, df_oxh_diff]
        titles  = ['OX GMI','OXM','OXH'] 
        # Create the plot
        if config.get('do_map_mae', False):
            plot_station_map_mae(df_list,spc,outdir=plot_outdir,titles=titles)
        if config.get('do_map_mb', False):
            plot_station_map_mb(df_list,spc,outdir=plot_outdir,titles=titles)

         
