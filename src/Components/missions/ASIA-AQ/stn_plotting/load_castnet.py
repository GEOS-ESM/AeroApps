#!/usr/bin/env python3
"""
    Script to load CASTNET observations and sampled model 
"""
import os, sys
import argparse
import yaml
import xarray as xr
import numpy as np
import pandas as pd
import argparse
import pytz

# map CASTNET timezone codes to standard timezone strings
tz_map = {
    'AL':   'US/Alaska',          # Alaska, with DST
    'AT':   'Canada/Atlantic',    # Atlantic, with DST
    'CE':   'US/Central',         # Central, with DST
    'EA':   'US/Eastern',         # Eastern, with DST
    'EAno': 'Etc/GMT+5',          # Eastern, no DST (fixed UTC-5)
    'HAno': 'Pacific/Honolulu',   # Hawaiian Aleutian, no DST
    'MO':   'US/Mountain',        # Mountain, with DST
    'MOno': 'Etc/GMT+7',          # Mountain, no DST (fixed UTC-7)
    'PA':   'US/Pacific',         # Pacific, with DST
}


def local_to_utc(row, col):
    """Convert a local time to UTC given a timezone string."""
    t = row[col]
    if pd.isna(t) or pd.isna(row['tz_str']):
        return pd.NaT
    try:
        tz       = pytz.timezone(row['tz_str'])
        t_local  = tz.localize(t, is_dst=None)   # is_dst=None raises error on ambiguous times
        t_utc    = t_local.astimezone(pytz.utc)
        return t_utc.replace(tzinfo=None)         # strip tzinfo for tz-naive output
    except pytz.exceptions.AmbiguousTimeError:
        # during DST fall-back there are ambiguous times -- assume DST=False
        tz      = pytz.timezone(row['tz_str'])
        t_local = tz.localize(t, is_dst=False)
        t_utc   = t_local.astimezone(pytz.utc)
        return t_utc.replace(tzinfo=None)
    except Exception as e:
        print(f"WARNING: Could not convert time {t} for site {row['SITE_ID']}: {e}")
        return pd.NaT



def get_mean_between_dates(ds, station, startDate, endDate, var='HNO3CONC'):
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

    # kg m-3 --> ug m-3
    ds_slice = ds_slice*1e9

    # Return the mean
    if len(ds_slice.time) == 0:
        return np.nan   # no data in this time window
    
    return float(ds_slice.mean(skipna=True).values)


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("config",help='configuration yaml file')

    args = parser.parse_args()
    model_name = 'ALL'

    config = yaml.safe_load(open(args.config))
    outdir = config['sampled_outdir']
    ctls = {
        'baseline': None,
        'oxm':      None,
        'oxh':      None,
    }
    # get model sampled data [HNO3CONC kg m-3]
    if model_name == 'ALL': 
        ctls['baseline']     = config.get('model_baseline')
        ctls['oxh']          = config.get('model_oxh')
        ctls['oxm']          = config.get('model_oxm')
    else:
        ctls[model_name]  = config.get('model_'+ model_name)


    datasets = {
        'baseline': None,
        'oxm':      None,
        'oxh':      None,
    }
    for key in ctls:
        if ctls[key] is not None:
            datasets[key] = xr.open_dataset(f'{outdir}/amon.{ctls[key]}.nc4')
 
    # get CASTNET site locations
    site_path = config['castnet_site_map']
    sites = pd.read_table(site_path, sep=',')

    castnet_path = config['castnet_data_dir']
    df_castnet   = pd.read_table(castnet_path, sep=',')

    # HNO3 ug m-3
    df_castnet = df_castnet.rename(columns={'NHNO3': 'HNO3'})
    df_castnet['startDate'] = pd.to_datetime(df_castnet['DATEON'])  # local time
    df_castnet['endDate']   = pd.to_datetime(df_castnet['DATEOFF'])


    # keep only rows in df_amon where SITEID appears in sites
    n_before = len(df_castnet['SITE_ID'].unique())

    df_castnet_save = df_castnet.copy()
    df_castnet = df_castnet[df_castnet['SITE_ID'].isin(sites['SITE_ID'])]

    n_after = len(df_castnet['SITE_ID'].unique())

    removed = n_before - n_after
    print(f"Sites before: {n_before}")
    print(f"Sites after:  {n_after}")
    print(f"Removed:      {removed}")


    # show which sites were removed
    removed_sites = set(df_castnet_save['SITE_ID'].unique()) - set(sites['SITE_ID'].unique())
    if removed_sites:
        print(f"Removed sites: {sorted(removed_sites)}")


    # convert local time to UTC

    # merge timezone info from sites into df_castnet
    df_castnet = df_castnet.merge(
        sites[['SITE_ID', 'TIME_ZONE']],
        on='SITE_ID',
        how='left'
    )

    # check for any unmapped timezone codes and remove them
    unknown_tz = set(df_castnet['TIME_ZONE'].unique()) - set(tz_map.keys())
    if unknown_tz:
        n_before = len(df_castnet)
        df_castnet = df_castnet[~df_castnet['TIME_ZONE'].isin(unknown_tz)]
        n_after = len(df_castnet)
        print(f"Removed {n_before - n_after} rows with unknown timezone codes: {unknown_tz}")
        print(f"Rows remaining: {n_after}")


    # map timezone codes to standard strings
    df_castnet['tz_str'] = df_castnet['TIME_ZONE'].map(tz_map)

    # convert startDate and endDate to UTC
    df_castnet['startDate_utc'] = df_castnet.apply(
        lambda row: local_to_utc(row, 'startDate'), axis=1)

    df_castnet['endDate_utc'] = df_castnet.apply(
        lambda row: local_to_utc(row, 'endDate'), axis=1)

    # drop helper column
    df_castnet = df_castnet.drop(columns=['tz_str'])

    # get time limits of the model dataset
    if model_name == 'ALL':
        model_name = 'baseline'
    ds_start = pd.Timestamp(datasets[model_name].time.values[0])
    ds_end   = pd.Timestamp(datasets[model_name].time.values[-1])

    # mask rows where observation window falls within model dataset time limits
    valid = (
        (df_castnet['startDate_utc'] >= ds_start) &
        (df_castnet['endDate_utc']   <= ds_end)
    )
    df_castnet = df_castnet[valid].copy()

    # get the model mean over the time of the observation
#    for name, ds in datasets.items():
#        if ds is not None:
#            df_castnet[name] = df_castnet.apply(
#                lambda row: get_mean_between_dates(
#                    ds        = ds,
#                    station   = row['SITE_ID'],
#                    startDate = row['startDate_utc'],
#                    endDate   = row['endDate_utc'],
#                    var       = 'HNO3CONC'
#                ),
#                axis=1
#            )

    # add a time column that is the mid-point time
    # normalize sets the hour to 00
    df_castnet['time'] = (df_castnet['startDate_utc'] + (df_castnet['endDate_utc'] - df_castnet['startDate_utc']) / 2).dt.normalize()

    # Quality Rating Code:
    # QA_CODE,0,Level 0 validated data,
    # QA_CODE,1,Level 1 validated data,
    # QA_CODE,1X,Level 1 data that are suspect based on manual screening.  Data should be omitted from any statistical analyses.,
    # QA_CODE,2,Level 2 validated data,
    # QA_CODE,3,Level 3 validated data,

    # keep all rows where QA_CODE == 3
    df_castnet = df_castnet[df_castnet['QA_CODE'] == 3].reset_index(drop=True)

    #Data Codes Contained in File:
    #COLUMN_NAME,CODE,DESCRIPTION,VALIDITY
    #%_F,!,X and U apply,Valid
    #%_F,#,Both U and L flags apply,Valid
    #%_F,$,L and U apply,Valid
    #%_F,&,Ambient concentration data measured/collected coincident with data flagged A,Valid
    #%_F,<,Missing < 25% of hourly sampling period,Valid
    #%_F,A,Anomalous value resulting from local anthropogenic activity,Valid
    #%_F,I,Invalid chemistry data and/or less than 75% valid flow data,Invalid
    #%_F,L,Less than 90% but greater than or equal to 75% valid flow data,Valid
    #%_F,M,Missing or completely invalid flow data,Invalid
    #%_F,N,Sample not analyzed,Invalid
    #%_F,NULL,Valid data,Valid
    #%_F,R,Re-run value used in calculation,Valid
    #%_F,S,Both L and R flags apply,Valid
    #%_F,U,Undetected - value listed is the reporting limit corrected by flow volume,Valid
    #%_F,X,Data suspect due to impregnation solution,Valid

    # invalid flag codes
    invalid_flags = ['I', 'M', 'N']

    # remove rows where NHNO3_F is invalid
    df_castnet = df_castnet[~df_castnet['NHNO3_F'].isin(invalid_flags)]


#    return df_castnet, sites
