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


    # show which sites are removed
    removed_sites = set(df_amon['SITEID'].unique()) - set(sites['siteId'].unique())
    if removed_sites:
        print(f"Removed sites: {sorted(removed_sites)}")

    # keep only rows in df_amon where SITEID appears in sites
    n_before = len(df_amon['SITEID'].unique())

    df_amon = df_amon[df_amon['SITEID'].isin(sites['siteId'])]

    n_after = len(df_amon['SITEID'].unique())

    removed = n_before - n_after
    print(f"Sites before: {n_before}")
    print(f"Sites after:  {n_after}")
    print(f"Removed:      {removed}")

    sys.exit()
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

    #                   lonmin, lonmax, latmin, latmax
    regions = {'NW': [  -130,   -100,     40,     50],
               'SW': [  -130,   -100,     30,     40],
               'SE': [  -100,    -70,     25,     38],
               'NE': [   -80,    -60,     38,     50],
               'MW': [  -100,    -80,     38,     50]}

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

    # =========================================================================
    # Time series plot -- monthly median +/- 1 std dev by region
    # =========================================================================

    # add month column to df_amon
    df_amon['month'] = pd.to_datetime(df_amon['startDate']).dt.month

    # color for each model
    model_colors = {
        'baseline': 'blue',
        'oxm':      'green',
        'oxh':      'red',
    }

    # one figure per region
    for region_name, (lonmin, lonmax, latmin, latmax) in regions.items():

        # subset df_amon to this region using lat/lon already merged in
        df_reg = df_amon[
            (df_amon['longitude'] >= lonmin) &
            (df_amon['longitude'] <= lonmax) &
            (df_amon['latitude']  >= latmin) &
            (df_amon['latitude']  <= latmax)
        ].copy()

        if len(df_reg) == 0:
            print(f"WARNING: No data for region {region_name} -- skipping")
            continue

        n_sites = df_reg['SITEID'].nunique()
        print(f"Region {region_name}: {n_sites} sites, {len(df_reg)} rows")

        fig, ax = plt.subplots(figsize=(12, 5))

        # --- plot AMoN observations ---
        amon_monthly = (df_reg.groupby('month')['NH3']
                        .agg(median='median', std='std')
                        .reset_index())

        ax.plot(amon_monthly['month'],
                amon_monthly['median'],
                color='black',
                linewidth=2,
                label='AMoN NH3',
                zorder=5)

        ax.fill_between(amon_monthly['month'],
                        amon_monthly['median'] - amon_monthly['std'],
                        amon_monthly['median'] + amon_monthly['std'],
                        color='black',
                        alpha=0.15)

        # --- plot each model ---
        for model in models:

            model_monthly = (df_reg.groupby('month')[model]
                             .agg(median='median', std='std')
                             .reset_index())

            # skip if no data for this model in this region
            if model_monthly['median'].isna().all():
                print(f"  WARNING: No {model} data in {region_name} -- skipping")
                continue

            ax.plot(model_monthly['month'],
                    model_monthly['median'],
                    linewidth=2,
                    label=model)

            ax.fill_between(model_monthly['month'],
                            model_monthly['median'] - model_monthly['std'],
                            model_monthly['median'] + model_monthly['std'],
                            alpha=0.15)

        # --- format axes ---
        ax.set_xticks(range(1, 13))
        ax.set_xticklabels(['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun',
                            'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])
        ax.set_xlabel('Month')
        ax.set_ylabel('NH3 (ug m-3)')
        ax.set_ylim(0, 5)
        ax.set_title(f'Monthly Median NH3 -- Region {region_name} '
                     f'({n_sites} sites)',
                     fontsize=13, fontweight='bold')
        ax.legend(fontsize=10)
        ax.grid(True, alpha=0.3, linestyle='--')

        # --- add region box info ---
        ax.text(0.01, 0.99,
                f'Lon: {lonmin} to {lonmax}\nLat: {latmin} to {latmax}',
                transform=ax.transAxes,
                verticalalignment='top',
                fontsize=9,
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

        plt.tight_layout()
        plt.savefig(f'plots/nh3_timeseries_{region_name}.png',
                    dpi=150, bbox_inches='tight')
        plt.show()
        plt.close()
