#!/usr/bin/env python3
"""
    Script to create plots comparing a GEOS simulations to each other sampled at IMPROVE sites 
    Reads in the outputs from sample/improve_sampler.py
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
from ts_plot_gas import ts_plot

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
    spcs = ['NH3','HNO3CONC']
    units = ['kg kg-1','kg m-3']

    # --- Molecular weights (g/mol) ---
    M_air = 28.97   # dry air

    M_gas = {'NH3': 17.03,
             'HNO3CONC': 63.01}


    # get output directory
    plot_outdir = config['plot_outdir'] + '/' + baseline_ctl
    os.makedirs(plot_outdir,exist_ok=True)

    baseline_ds = xr.open_dataset(f'{outdir}/improve.{baseline_ctl}.nc4')
    oxh_ds      = xr.open_dataset(f'{outdir}/improve.{oxh_ctl}.nc4')
    oxm_ds      = xr.open_dataset(f'{outdir}/improve.{oxm_ctl}.nc4')

    # get IMPROVE site locations
    site_path = config['improve_site_map']
    sites = SITE_MAP(site_path)
                    # lonmin, lonmax, latmin, latmax
    regions = {'NW': [  -130,   -100,     40,     50],
               'SW': [  -130,   -100,     30,     40],
               'SE': [  -100,    -70,     25,     38],
               'NE': [   -80,    -60,     38,     50],
               'MW': [  -100,    -80,     38,     50]}
    for reg in regions:
        lon_min, lon_max, lat_min, lat_max = regions[reg]

        for spc in spcs:
            # Get baseline data
            ds_baseline = baseline_ds[spc]
            ds_oxh = oxh_ds[spc]
            ds_oxm = oxm_ds[spc]

            # convert to pptv
            if spc == 'NH3':
                ds_baseline = ds_baseline*baseline_ds['AIRDENS']
                ds_ohx      = ds_oxh*oxh_ds['AIRDENS']
                ds_oxm      = ds_oxm*oxm_ds['AIRDENS']

            mw_gas = M_gas[spc]
            # kg/m3 -> mol/m3 (divide by molecular weight in kg/mol)
            ds_baseline = ds_baseline / (mw_gas / 1000.0)        
            ds_oxh = ds_oxh / (mw_gas / 1000.0)
            ds_oxm = ds_oxm / (mw_gas / 1000.0)

            # mol/m3 -> mol/mol 
            # air molar density = airdens / (M_air/1000) mol/m3
            air_mol_m3 = baseline_ds['AIRDENS'] / (M_air / 1000.0)
            ds_baseline = ds_baseline / air_mol_m3

            air_mol_m3 = oxh_ds['AIRDENS'] / (M_air / 1000.0)
            ds_oxh = ds_oxh / air_mol_m3

            air_mol_m3 = oxm_ds['AIRDENS'] / (M_air / 1000.0)
            ds_oxm = ds_oxm / air_mol_m3        
           
            # mol/mol -> pptv
            ds_baseline = ds_baseline * 1e12
            ds_oxh = ds_oxh * 1e12
            ds_oxm = ds_oxm * 1e12
             
            # Create dataframes for baseline vs obs (top subplot)
            hourly_mean = ds_baseline.sel(lev=72).groupby('time.hour').mean(skipna=True).rename({'hour': 'time'})
            df_baseline = hourly_mean.to_dataframe(name=spc).reset_index()
            df_baseline["Source"] = "Baseline"
            # Merge lat/lon from sites.df into df_baseline_diff
            # matching on SiteCode (sites.df) = station (df_baseline_diff)
            df_baseline = df_baseline.merge(
                sites.df[['SiteCode', 'Latitude', 'Longitude']],
                left_on='station',
                right_on='SiteCode',
                how='left'
            ).drop(columns='SiteCode')  # drop the redundant SiteCode column
            mask = ((df_baseline['Latitude']  >= lat_min) & (df_baseline['Latitude']  <= lat_max) &
                    (df_baseline['Longitude'] >= lon_min) & (df_baseline['Longitude'] <= lon_max))

            df_baseline = df_baseline[mask]
            
            hourly_mean = ds_oxh.sel(lev=72).groupby('time.hour').mean(skipna=True).rename({'hour': 'time'})
            df_oxh = hourly_mean.to_dataframe(name=spc).reset_index()
            df_oxh["Source"] = 'OxH'
            # Merge lat/lon from sites.df into df_baseline_diff
            # matching on SiteCode (sites.df) = station (df_baseline_diff)
            df_oxh = df_oxh.merge(
                sites.df[['SiteCode', 'Latitude', 'Longitude']],
                left_on='station',
                right_on='SiteCode',
                how='left'
            ).drop(columns='SiteCode')  # drop the redundant SiteCode column
            mask = ((df_oxh['Latitude']  >= lat_min) & (df_oxh['Latitude']  <= lat_max) &
                    (df_oxh['Longitude'] >= lon_min) & (df_oxh['Longitude'] <= lon_max))

            df_oxh = df_oxh[mask]


            hourly_mean = ds_oxm.sel(lev=72).groupby('time.hour').mean(skipna=True).rename({'hour': 'time'})
            df_oxm = hourly_mean.to_dataframe(name=spc).reset_index()
            df_oxm["Source"] = 'OxM'
            # Merge lat/lon from sites.df into df_baseline_diff
            # matching on SiteCode (sites.df) = station (df_baseline_diff)
            df_oxm = df_oxm.merge(
                sites.df[['SiteCode', 'Latitude', 'Longitude']],
                left_on='station',
                right_on='SiteCode',
                how='left'
            ).drop(columns='SiteCode')  # drop the redundant SiteCode column
            mask = ((df_oxm['Latitude']  >= lat_min) & (df_oxm['Latitude']  <= lat_max) &
                    (df_oxm['Longitude'] >= lon_min) & (df_oxm['Longitude'] <= lon_max))

            df_oxm = df_oxm[mask]

            
            # Create dataframe for top subplot (baseline vs obs)
            plot_df_top = pd.concat([df_baseline, df_oxh, df_oxm], ignore_index=True)
            
            # Create the plot
            unit = 'pptv'
            os.makedirs(plot_outdir+'/'+reg+'_hourly', exist_ok=True)
            ts_plot(plot_df_top, spc, unit, plot_outdir+'/'+reg+'_hourly')

             
