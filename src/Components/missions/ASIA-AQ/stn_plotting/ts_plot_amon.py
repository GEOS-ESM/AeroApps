#!/usr/bin/env python3
"""
    Script to create plots comparing a GEOS simulations to each other sampled at AMoN sites 
    Reads in the outputs from sample/amon_sampler.py
"""
import os, sys
import argparse
import yaml
import xarray as xr
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats
from map_plot import plot_station_map

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
    spcs = ['NH3']
    units = ['ug m-3']

    # --- Molecular weights (g/mol) ---
    M_air = 28.97   # dry air

    M_gas = 17.03


    # get output directory
    plot_outdir = config['plot_outdir'] + '/' + baseline_ctl
    os.makedirs(plot_outdir,exist_ok=True)

    baseline_ds = xr.open_dataset(f'{outdir}/amon.{baseline_ctl}.nc4')
    oxh_ds      = xr.open_dataset(f'{outdir}/amon.{oxh_ctl}.nc4')
    oxm_ds      = xr.open_dataset(f'{outdir}/amon.{oxm_ctl}.nc4')

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

    df_amon = df_amon[df_amon['SITEID'].isin(sites['siteId'])]

    n_after = len(df_amon['SITEID'].unique())

    removed = n_before - n_after
    print(f"Sites before: {n_before}")
    print(f"Sites after:  {n_after}")
    print(f"Removed:      {removed}")

    # show which sites were removed
    removed_sites = set(df_amon['SITEID'].unique()) - set(sites['siteId'].unique())
    if removed_sites:
        print(f"Removed sites: {sorted(removed_sites)}")


    datasets = {
        'Baseline': baseline_ds,
        'OXM':      oxm_ds,
        'OXH':      oxh_ds,
    }

    # filter out values greater than 3 STD
    # rationale is these represent local sources that
    # aren't represented in the model
    mean = df_amon['NH3'].mean()
    std  = df_amon['NH3'].std()

    threshold_3sigma = mean + 3 * std
    valid = df_amon['NH3'] < threshold_3sigma
    df_amon = df_amon[valid].copy()


    # get time limits of this dataset
    ds_start = pd.Timestamp(baseline_ds.time.values[0])
    ds_end   = pd.Timestamp(baseline_ds.time.values[-1])

    # mask rows where observation window falls within dataset time limits
    valid = (
        (df_amon['startDate'] >= ds_start) &
        (df_amon['endDate']   <= ds_end)
    )
    df_amon = df_amon[valid].copy()

    # get the model mean over the time of the observation
    for name, ds in datasets.items():

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

    # merge lat/lon from sites into df_amon on site ID
    df_amon = df_amon.merge(
        sites[['siteId', 'latitude', 'longitude']],
        left_on='SITEID',
        right_on='siteId',
        how='left'
    )


                    # lonmin, lonmax, latmin, latmax
    regions = {'NW': [  -130,   -100,     40,     50],
               'SW': [  -130,   -100,     30,     40],
               'SE': [  -100,    -70,     25,     38],
               'NE': [   -80,    -60,     38,     50],
               'MW': [  -100,    -80,     38,     50]}

    region_name = 'SE'
    lonmin, lonmax, latmin, latmax = regions[region_name]

    # apply spatial filter
    df_region = df_amon[
        (df_amon['longitude'] >= lonmin) &
        (df_amon['longitude'] <= lonmax) &
        (df_amon['latitude']  >= latmin) &
        (df_amon['latitude']  <= latmax)
    ]



    # make 2D KDE plots
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))

    # get common axis limits across all three panels
    models = list(datasets.keys())
    df_plot = df_region[['NH3'] + models].dropna(how='all')
    xmax = df_plot['NH3'].quantile(0.99)
    ymax = df_plot[models].quantile(0.99).max()
    lim  = max(xmax, ymax)

    for ax, model in zip(axes,models):

        df_sub = df_region[['NH3', model]].dropna()
        diff   = df_sub[model] - df_sub['NH3']

        # 1D KDE of difference
        diff.plot.kde(ax=ax, linewidth=2, label='KDE')

        # reference lines
        ax.axvline(0,           color='r', linestyle='--', linewidth=1.5,
                   label='zero')
        ax.axvline(diff.mean(), color='k', linestyle='--', linewidth=1.5,
                   label=f'mean={diff.mean():.2f}')

        # stats
        n    = len(diff)
        std  = diff.std()
        rmse = np.sqrt((diff**2).mean())

        ax.text(0.05, 0.95,
                f'Bias={diff.mean():.2f}\nStd={std:.2f}\nRMSE={rmse:.2f}\nN={n}',
                transform=ax.transAxes,
                verticalalignment='top',
                fontsize=10,
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

        ax.set_xlabel('Model - AMoN NH3 (ug/m3)')
        ax.set_ylabel('Density')
        ax.set_title(model)
        ax.legend()

    plt.suptitle('Model minus AMoN NH3', fontsize=14, y=1.02)
    plt.tight_layout()
#    plt.savefig('nh3_diff_kde.png', dpi=150, bbox_inches='tight')
    plt.show()
