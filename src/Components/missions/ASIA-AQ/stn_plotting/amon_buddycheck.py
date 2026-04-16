#!/usr/bin/env python3
"""
    Script to run a  buddy check on AMoN sites 
"""
import os, sys
import argparse
import yaml
import xarray as xr
import numpy as np
import pandas as pd
from pyobs import buddycheck

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
    df_amon['model'] = df_amon.apply(
        lambda row: get_mean_between_dates(
            ds        = baseline_ds,
            station   = row['SITEID'],
            startDate = row['startDate'],
            endDate   = row['endDate'],
            var       = 'NH3'
        ),
        axis=1
    )

    # rename columns to names expected by buddycheck
    df_amon = df_amon.rename(columns={
            'SITEID': 'station_id',
            'NH3': 'obs',
            })

    sites = sites.rename(columns={
            'siteId': 'station_id',
            'latitude': 'lat',
            'longitude': 'lon',
            })

    # add a time column that is the mid-point time
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

    # configuration parameters for the buddy check
    config = buddycheck.InnovationBuddyConfig(
                radius_km=150.0,
                min_neighbors=3,
                distance_scale_km=75.0,
                buddy_method="weighted_mean",
                outlier_sigma_thresh=2.5,
                min_valid_buddy_count=10,
    )

    buddy_time_df, station_scores_df = buddycheck.score_stations_innovation_based(
        obs_df=df_amon,
        station_meta=sites,
        config=config,
    )

    # implement an alternative review flag
    review_flag = ((station_scores_df["suspect_score"] >= 6.0) &
                  (
                    (station_scores_df["z_f_bias_vs_buddies"] >= 2.5) |
                    (station_scores_df["z_f_rmse_vs_buddies"] >= 2.5) |
                    (station_scores_df["z_f_frac_outlier"] >= 2.5)
                  )
                 )

    station_scores_df["review_flag_new"] = review_flag

    # print the top 20 stations
    print(station_scores_df.sort_values("suspect_score", ascending=False).head(20))

    files = buddycheck.plot_flagged_stations(
        station_scores_df=station_scores_df,
        buddy_time_df=buddy_time_df,
        station_meta=sites,
        outdir="plots/amon_flagged_stations",
        only_flagged=False,
        top_n=5,
        )
