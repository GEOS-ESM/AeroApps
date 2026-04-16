#!/usr/bin/env python3
"""
    Script to load AMoN observations and sampled model 
"""
import os, sys
import argparse
import yaml
import xarray as xr
import numpy as np
import pandas as pd

def get_mean_between_dates(ds, station, startDate, endDate, var='NH3'):
    """
    Get the mean of a variable in an xarray Dataset 
    between startDate and endDate for a given station.

    Parameters
    ----------
    ds        : xarray Dataset
    station   : station ID string
    startDate : start date (datetime)
    endDate   : end date (datetime)
    var       : variable name to average

    Returns
    -------
    mean value between the two dates
    """
    # Select the station and time window
    ds_slice = ds[var].sel(
        lev=72,
        station=station,
        time=slice(startDate, endDate)
    )

    ds_airdens = ds['AIRDENS'].sel(
        lev=72,
        station=station,
        time=slice(startDate, endDate)
    )

    # kg kg-1 --> ug m-3
    ds_slice = ds_slice*ds_airdens*1e9

    # Return the mean
    if len(ds_slice.time) == 0:
        return np.nan   # no data in this time window
    
    return float(ds_slice.mean(skipna=True).values)


def LOAD_AMON(config,model_name='ALL'):

    config = yaml.safe_load(open(config))

    ctls = {
        'baseline': None,
        'oxm':      None,
        'oxh':      None,
    }
    # get model sampled data
    if model_name == 'ALL': 
        ctls['baseline']     = config.get('model_baseline')
        ctls['oxh']          = config.get('model_oxh')
        ctls['oxm']          = config.get('model_oxm')
    else:
        ctls[model_name]  = config.get('model_'+ model_name)


    datasets = {
        'baseline': baseline_ds,
        'oxm':      oxm_ds,
        'oxh':      oxh_ds,
    }
    for key in ctls:
        if ctls[key] is not None:
            datasets[key] = xr.open_dataset(f'{outdir}/amon.{ctls[key]}.nc4')
 
    # get AMoN site locations
    site_path = config['amon_site_map']
    sites = pd.read_table(site_path, sep=',')

    amon_path = config['amon_data_dir']
    df_amon   = pd.read_table(amon_path, sep=',')
    df_amon = df_amon.rename(columns={'CONC': 'NH3'})
    df_amon['startDate'] = pd.to_datetime(df_amon['startDate'])
    df_amon['endDate']   = pd.to_datetime(df_amon['endDate'])


    # keep only rows in df_amon where SITEID appears in sites
    n_before = len(df_amon['SITEID'].unique())

    df_amon_save = df_amon.copy()
    df_amon = df_amon[df_amon['SITEID'].isin(sites['siteId'])]

    n_after = len(df_amon['SITEID'].unique())

    removed = n_before - n_after
    print(f"Sites before: {n_before}")
    print(f"Sites after:  {n_after}")
    print(f"Removed:      {removed}")


    # show which sites were removed
    removed_sites = set(df_amon_save['SITEID'].unique()) - set(sites['siteId'].unique())
    if removed_sites:
        print(f"Removed sites: {sorted(removed_sites)}")


    # get time limits of the model dataset
    if model_name == 'ALL':
        model_name = 'baseline'
    ds_start = pd.Timestamp(datasets[model_name].time.values[0])
    ds_end   = pd.Timestamp(datasets[model_name].time.values[-1])

    # mask rows where observation window falls within model dataset time limits
    valid = (
        (df_amon['startDate'] >= ds_start) &
        (df_amon['endDate']   <= ds_end)
    )
    df_amon = df_amon[valid].copy()

    # get the model mean over the time of the observation
    for name, ds in datasets.items():
        if ds is not None:
            df_amon[name] = df_amon.apply(
                lambda row: get_mean_between_dates(
                    ds        = ds,
                    station   = row['SITEID'],
                    startDate = row['startDate'],
                    endDate   = row['endDate'],
                    var       = 'NH3'
                ),
                axis=1
            )

    # add a time column that is the mid-point time
    # normalize sets the hour to 00
    df_amon['time'] = (df_amon['startDate'] + (df_amon['endDate'] - df_amon['startDate']) / 2).dt.normalize()

    # Quality Rating Code:
    # A: valid data
    # B: valid data with minor problems
    # C: invalid data
    # keep all rows where QR != C
    df_amon = df_amon[df_amon['QR'] != 'C'].reset_index(drop=True)

    # Replicate
    # Sample replicate (A, B, C) or travel blank (T)
    # keep all rows where REPLICATE != T
    df_amon = df_amon[df_amon['REPLICATE'] != 'T'].reset_index(drop=True)


    return df_amon
