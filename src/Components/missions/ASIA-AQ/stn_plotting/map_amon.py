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
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import cartopy.crs as ccrs
import cartopy.feature as cfeature


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


    # get station means
    models = list(datasets.keys())
    station_means = []
    for model in models:

        df_sub = df_amon[['NH3', model]].dropna()
        diff   = df_sub[model] - df_sub['NH3']
        df_amon['diff_'+model] = np.nan
        df_amon.loc[df_sub.index, 'diff_'+model] = diff

    
        stn_mean = (df_amon.groupby(['SITEID', 'latitude', 'longitude'])['diff_'+model]
                      .apply(lambda x: np.nanmean(x))
                      .reset_index()
                      .rename(columns={'diff_'+model: 'abs_mean'}))

        station_means.append(stn_mean)


    # --- Set colorbar limits ---
    vmin = -5
    vmax = 5

    colors = ['#08306B',   # dark blue
         '#2171B5',   # medium blue
         '#6BAED6',   # light blue
         '#74C476',   # light green
         '#238B45',   # medium green
         '#FFFFFF',   # white (center)
         '#9E9AC8',   # light purple
         '#6A51A3',   # medium purple
         '#F16913',   # light red/orange
         '#D62728',   # medium red
         '#67000D']  # dark red
    cmap = mcolors.LinearSegmentedColormap.from_list(
        'blue_green_white_purple_red', colors, N=256)


    # --- Compute map extent across all 3 dataframes ---
    all_lons = np.concatenate([sm['longitude'].values for sm in station_means])
    all_lats = np.concatenate([sm['latitude'].values for sm in station_means])
    lon_pad  = 5
    lat_pad  = 5
    extent   = [all_lons.min() - lon_pad,
                #all_lons.max() + lon_pad,
                -60,
                all_lats.min() - lat_pad,
                all_lats.max() + lat_pad]



    # make 2D KDE plots
    fig, axes = plt.subplots(3, 1, figsize=(9, 12),subplot_kw={'projection': ccrs.PlateCarree()})


    for idx, (ax, station_mean, model) in enumerate(zip(axes, station_means, models)):

        # --- Add map features ---
        ax.add_feature(cfeature.LAND,       facecolor='lightgray', alpha=0.3)
        ax.add_feature(cfeature.OCEAN,      facecolor='lightblue', alpha=0.3)
        ax.add_feature(cfeature.COASTLINE,  linewidth=0.5)
        ax.add_feature(cfeature.BORDERS,    linewidth=0.5, linestyle='--')
        ax.add_feature(cfeature.STATES,     linewidth=0.3, linestyle=':')

        # --- Set map extent ---
        ax.set_extent(extent, crs=ccrs.PlateCarree())

        # --- Plot scatter points ---
        sc = ax.scatter(station_mean['longitude'],
                        station_mean['latitude'],
                        c=station_mean['abs_mean'],
                        cmap=cmap,
                        vmin=vmin,
                        vmax=vmax,
                        s=40,
                        edgecolors='black',
                        linewidths=0.5,
                        transform=ccrs.PlateCarree(),
                        zorder=5)

        # --- Add gridlines ---
        gl = ax.gridlines(draw_labels=True, linewidth=0.5,
                          color='gray', alpha=0.5, linestyle='--')
        gl.top_labels   = False
        gl.right_labels = False

        # --- Add title ---
        ax.set_title(model, fontsize=13, fontweight='bold', pad=10)


    # --- After the loop, add colorbar to the right of the last axis ---
    # Create a new axes on the right side of the figure
    fig.subplots_adjust(right=0.85)  # make room on the right

    # Add a new axes for the colorbar
    cbar_ax = fig.add_axes([0.88,   # left position
                             0.1,   # bottom position
                             0.02,  # width
                             0.8])  # height

    cbar = fig.colorbar(sc, cax=cbar_ax, orientation='vertical')
    cbar.set_label(f'Mean Bias NH3 (ug m-3)', fontsize=11)


#    plt.suptitle('Model minus AMoN NH3', fontsize=14, y=1.02)
#    plt.tight_layout()
#    plt.savefig('nh3_diff_kde.png', dpi=150, bbox_inches='tight')
    plt.show()
