#!/usr/bin/env python3

'''
This code subsamples a single GEOS model experiment according to hourly means of AERONET observations.
The output is a CSV file with high-level statistics for all sites with observations for the 
requested date range including means, standard deviations, RMSE, and correlations as well as
a time series CSV file at the requested temporal frequency (hourly, daily, monthly).
Since this utilizes aeronet.py from GMAOpyobs, the python path must be se first (export PYTHONPATH="/discover/nobackup/acollow/GMAOpyobs/src/:$PYTHONPATH")
'''

import os
import numpy as np
import pandas as pd
import xarray as xr
from datetime import datetime, timedelta
import warnings
import multiprocessing as mp
from functools import partial
import time
from pyobs.aeronet import AERONET_L2, VARS 
import sys
import yaml

warnings.filterwarnings('ignore')


#####################################################################################
# Parse the aeronet daily files
#####################################################################################
def parse_aeronet_daily_file(filepath, start_date, end_date):
    """Reads a daily global AERONET file using AERONET_L2 class, applies QC, and returns hourly means per station."""
    try:
        custom_vars = list(VARS) + ['440_870_Angstrom_Exponent']
         
        aero = AERONET_L2(filepath, version=3, Vars=custom_vars, Verbose=False)
     
        if aero.nobs == 0:
            return None, "No data loaded"
        
        df = pd.DataFrame({
            'DateTime': aero.tyme,
            'lats': aero.Latitude,
            'lons': aero.Longitude,
            'station': aero.Location, 
            'aod_550': aero.AOT_550,
            'angstrom': getattr(aero, '440_870_Angstrom_Exponent', np.nan) 
        })

        df = df[(df['DateTime'] >= start_date) & (df['DateTime'] <= end_date)]
        if df.empty: 
            return None, "No data in specified date range"
        
        df.replace(-999.0, np.nan, inplace=True)

        df.loc[(df['aod_550'] < 0) | (df['aod_550'] >= 10), 'aod_550'] = np.nan

        valid_mask_870 = (df['angstrom'] >= -1) & (df['angstrom'] <= 3)
        df.loc[~valid_mask_870, 'angstrom'] = np.nan

        df = df.dropna(subset=['aod_550', 'angstrom'])

        if df.empty:
            return None, "No valid paired AOD/Angstrom data after QC"

        df['hour'] = df['DateTime'].dt.floor('H')
        
        hourly = df.groupby(['station', 'hour']).agg({
            'aod_550': 'mean', 
            'angstrom': 'mean',
            'lats': 'first',
            'lons': 'first'
        }).reset_index()
        
        return hourly, "Success"
        
    except ValueError as e:
        return None, f"Parsing Error: {str(e)}"
    except Exception as e:
        return None, f"Exception test: {type(e).__name__} - {str(e)}"
#####################################################################################
# Process the analysis data
#####################################################################################

def processanalysis(date_tuple, analysis_path_template):

    date, daily_df = date_tuple

    for hour_val, group in daily_df.groupby('hour'):
        YYYYMM = hour_val.strftime('%Y%m')
        YYYYMMDD = hour_val.strftime('%Y%m%d')
        MM = hour_val.strftime('%m')
        
        HH = hour_val.strftime('%H')
       
        analysis_file = analysis_path_template.format(YYYYMM=YYYYMM, YYYYMMDD=YYYYMMDD, HH=HH, MM=MM)
        if os.path.exists(analysis_file):
            with xr.open_dataset(analysis_file) as ds_mod:
                for idx, row in group.iterrows():
                    mlon = row.lons + 360 if (ds_mod.lon.min() >= 0 and row.lons < 0) else row.lons
                    pt = ds_mod.sel(lat=row['lats'], lon=mlon, method='nearest')
                    if 'time' in pt.dims:
                        pt = pt.isel(time=0)
                    daily_df.at[idx, 'analysis_aod'] = float(pt['TOTEXTTAU'].values) # Changed from 'TOTEXTTAU' to 'BCEXTTAU' because 'TOTEXTTAU' was not in dataset
                    daily_df.at[idx, 'analysis_ang'] = float(pt['TOTANGSTR'].values) # Changed from 'TOTANGSTR' to 'NIEXTTAU' because 'TOTANGSTR' was not in dataset 

    return daily_df


#####################################################################################
# Process the model data
#####################################################################################
def processmodel(date_tuple, model_path_template):
    """Processes all stations for a single day against the GEOS model data"""
    date, daily_df = date_tuple

    for hour_val, group in daily_df.groupby('hour'):
        YYYYMM = hour_val.strftime('%Y%m')
        YYYYMMDD = hour_val.strftime('%Y%m%d')
        
        HH = hour_val.strftime('%H')
       
        model_file = model_path_template.format(YYYYMM=YYYYMM, YYYYMMDD=YYYYMMDD, HH=HH)
        if os.path.exists(model_file):
            with xr.open_dataset(model_file) as ds_mod:
                for idx, row in group.iterrows():

                    mlon = row['lons'] + 360 if (ds_mod.lons.min() >= 0 and row['lons'] < 0) else row['lons']
                    pt = ds_mod.sel(Ydim=row['lats'], Xdim=mlon, method='nearest')
                    if 'time' in pt.dims:
                        pt = pt.isel(time=0)
                    daily_df.at[idx, 'model_aod'] = float(pt['TOTEXTTAU'].values) # Changed from 'TOTEXTTAU' to 'BCEXTTAU' because 'TOTEXTTAU' was not in dataset
                    daily_df.at[idx, 'model_ang'] = float(pt['TOTANGSTR'].values) # Changed from 'TOTANGSTR' to 'NIEXTTAU' because 'TOTANGSTR' was not in dataset 

    return daily_df


#####################################################################################
# Main function
#####################################################################################
def main(start_date, end_date, base_path, experiment_name, output_dir="./aeronet_model_comparison/", min_points=50, ts_freq='none'):
    start_time = time.time()
    n_procs = max(1, mp.cpu_count() - 1)
    #os.makedirs(output_dir, exist_ok=True)
    model_path_template = os.path.join(config['paths']['forecast'], f"CYCLED_REPLAY_P10800_C21600_T21600_{{YYYYMMDD}}_00z", f"GEOS.hwt_15mn_slv_LCC.{{YYYYMMDD}}_{{HH}}00z.nc4")
    
    aeronet_dir_base = config['paths']['aeronet']

    analysis_path_template = os.path.join(config['paths']['fp_forecast'], f'M{{MM}}', f'f5430_fp.inst1_2d_hwl_Nx.{{YYYYMMDD}}_{{HH}}00z.nc4')

    file_comparison = 'aeronet'

    print(f" Locating AERONET daily files between {start_date.strftime('%Y-%m-%d')} and {end_date.strftime('%Y-%m-%d')}...")
    
    files = []
    current_date = start_date.replace(hour=0, minute=0, second=0)


    ### Processing aeronet files


    while current_date <= end_date:
        YYYY = current_date.strftime('%Y')
        MM = current_date.strftime('%m')
        YYYYMMDD = current_date.strftime('%Y%m%d')
        #Aeronet switches filename formats in 2024!
        possible_filenames = [
            f"aeronet_v30.{YYYYMMDD}.txt",
            f"aeronet_v3.aod.{YYYYMMDD}.txt"
        ]
        
        for fname in possible_filenames:
            file_path = f"{aeronet_dir_base}/Y{YYYY}/M{MM}/{fname}"
            
            if os.path.exists(file_path):
                files.append(file_path)
                break 
        
        current_date += timedelta(days=1)
    
    if not files:
        print(f"No AERONET files found in the specified date range at {aeronet_dir_base}!")
        return
        
    print(f"Found {len(files)} daily AERONET files. Parsing...")
        
    parse_func = partial(parse_aeronet_daily_file, start_date=start_date, end_date=end_date)
    dfs = []
    error_counts = {}

    with mp.Pool(n_procs) as pool:
        for df, msg in pool.imap_unordered(parse_func, files):
            if df is not None:
                dfs.append(df)
            else:
                print(msg)
                error_counts[msg] = error_counts.get(msg, 0) + 1

    print("\n--- AERONET Parse Summary ---")
    print(f"Successfully loaded {len(dfs)} days of data.")
    
    if error_counts:
        print("\nReasons for skipped files:")
        for reason, count in sorted(error_counts.items(), key=lambda x: x[1], reverse=True):
            print(f"  - {reason}: {count} files")
            
    if not dfs:
        print("\nNo data to process. Exiting.")
        return
        
    master_df = pd.concat(dfs, ignore_index=True)
    master_df['date'] = master_df['hour'].dt.floor('D')
    print(f"\nTotal valid AERONET hours across all stations: {len(master_df)}")
    
    print(f"\nSampling {experiment_name}...")
    master_df['model_aod'] = np.nan
    master_df['model_ang'] = np.nan
    master_df['analysis_aod'] = np.nan
    master_df['analysis_ang'] = np.nan

    grouped = [(date, group) for date, group in master_df.groupby('date')]


    ### Processing model data


    if "forecast" in config['evaluation']:
        file_comparison = file_comparison + '_fcst'
    #grouped = [(date, group) for date, group in master_df.groupby('date')]
        process_func = partial(processmodel, model_path_template=model_path_template)
        processed_dfs = []
    
        with mp.Pool(n_procs) as pool:
            for i, result_df in enumerate(pool.imap_unordered(process_func, grouped)):
                processed_dfs.append(result_df)
                if (i+1) % 10 == 0 or (i+1) == len(grouped):
                    print(f"Processed {i+1}/{len(grouped)} model days...")
        model_df = pd.concat(processed_dfs, ignore_index=True)
    
        final_df = model_df.dropna(subset=['model_aod', 'model_ang'])


    ### Processing analysis data


    if "analysis" in config['evaluation']:
        
        file_comparison = file_comparison + '_analysis'
        analysis_process_func = partial(processanalysis, analysis_path_template=analysis_path_template)
        analysis_processed_dfs = []


        with mp.Pool(n_procs) as pool:
            for i, result_df in enumerate(pool.imap_unordered(analysis_process_func, grouped)):
                analysis_processed_dfs.append(result_df)
                if (i+1) % 10 == 0 or (i+1) == len(grouped):
                    print(f"Processed {i+1}/{len(grouped)} analysis days...")

        analysis_df = pd.concat(analysis_processed_dfs, ignore_index = True)
        analysis_df = analysis_df.dropna(subset=['analysis_aod', 'analysis_ang'])
        try:
            final_df = model_df.combine_first(analysis_df)
        except:
            final_df = analysis_df
    valid_stations = final_df['station'].value_counts()[final_df['station'].value_counts() >= min_points].index

    final_df = final_df[final_df['station'].isin(valid_stations)].copy()


    ### Calculating Stats


    print(f"\nCalculating Statistics... ({len(final_df)} matched hour-station points from {len(valid_stations)} stations)")
    stats = []
    for station, grp in final_df.groupby('station'):
        stats.append({
            'Experiment': experiment_name,
            'Station': station,
            'Lats': grp['lats'].iloc[0],
            'Lons': grp['lons'].iloc[0],
            'N_Points': len(grp),
            
            'AERONET_AOD_Mean': np.nanmean(grp['aod_550']),
            'AERONET_AOD_Std': np.nanstd(grp['aod_550']),
            'MODEL_AOD_Mean': np.nanmean(grp['model_aod']),
            'MODEL_AOD_Std': np.nanstd(grp['model_aod']),
            'ANALYSIS_AOD_Mean': np.nanmean(grp['analysis_aod']),
            'ANALYSIS_AOD_Std': np.nanstd(grp['analysis_aod']),

            
            'AERONET_Angstrom_Mean': np.nanmean(grp['angstrom']),
            'AERONET_Angstrom_Std': np.nanstd(grp['angstrom']),
            'MODEL_Angstrom_Mean': np.nanmean(grp['model_ang']),
            'MODEL_Angstrom_Std': np.nanstd(grp['model_ang']),
            'ANALYSIS_Angstrom_Mean': np.nanmean(grp['analysis_ang']),
            'ANALYSIS_Angstrom_Std': np.nanstd(grp['analysis_ang']),
            
            'MODEL_AOD_RMSE': np.sqrt(((grp['model_aod'] - grp['aod_550'])**2).mean()),
            'MODEL_Angstrom_RMSE': np.sqrt(((grp['model_ang'] - grp['angstrom'])**2).mean()),
            'ANALYSIS_AOD_RMSE': np.sqrt(((grp['analysis_aod'] - grp['aod_550'])**2).mean()),
            'ANALYSIS_Angstrom_RMSE': np.sqrt(((grp['analysis_ang'] - grp['angstrom'])**2).mean()),
            
            'MODEL_AOD_Correlation': grp['model_aod'].corr(grp['aod_550']),
            'MODEL_Angstrom_Correlation': grp['model_ang'].corr(grp['angstrom']),
            'ANALYSIS_AOD_Correlation': grp['analysis_aod'].corr(grp['aod_550']),
            'ANALYSIS_Angstrom_Correlation': grp['analysis_ang'].corr(grp['angstrom'])

        })
       



    date_str = f"{start_date.strftime('%Y%m%d')}_{end_date.strftime('%Y%m%d')}"
    
    if stats:
        out_df = pd.DataFrame(stats)  
        output_file = os.path.join(output_dir[:19],output_dir[28:], f"{experiment_name}_{file_comparison}_comparison_stats_{ts_freq}_{date_str}.csv")
        out_df.to_csv(output_file, index=False)
        print(f"Stats saved to {output_file}")
    else:
        print(f"No stations met the minimum threshold of {min_points} valid matched points to calculate statistics.")

    if ts_freq != 'none' and not final_df.empty:
        print(f"\nGenerating {ts_freq} timeseries data...")
        
        cols_to_keep = ['hour', 'station', 'lats', 'lons', 'aod_550', 'model_aod', 'analysis_aod', 'angstrom', 'model_ang', 'analysis_ang']
        ts_df = final_df[cols_to_keep].copy()
        ts_df.rename(columns={'hour': 'time', 'aod_550': 'aeronet_aod', 'angstrom': 'aeronet_angstrom'}, inplace=True)
        
        ts_df.insert(0, 'experiment', experiment_name)
        
        if ts_freq == 'daily':
            ts_df['time'] = ts_df['time'].dt.floor('D')
            ts_df = ts_df.groupby(['experiment', 'station', 'time', 'lats', 'lons']).mean(numeric_only=True).reset_index()
        elif ts_freq == 'monthly':
            ts_df['time'] = ts_df['time'].dt.to_period('M').dt.to_timestamp()
            ts_df = ts_df.groupby(['experiment', 'station', 'time', 'lats', 'lons']).mean(numeric_only=True).reset_index()
            
        ts_output_file = os.path.join(output_dir[:19],output_dir[28:], f"{experiment_name}_{file_comparison}_timeseries_{ts_freq}_{date_str}.csv")
        
        ts_df.to_csv(ts_output_file, index=False)
        print(f"Time series data saved to {ts_output_file}")

    print(f"Total time: {(time.time() - start_time)/60:.1f} minutes!")


#####################################################################################
# Running main
#####################################################################################
if __name__ == "__main__":
    # to change settings for the date, temporal frequency, and what is beign evaluatated, edit config.yaml

    with open('config.yaml', 'r') as file:
        config = yaml.safe_load(file)
    start_date = datetime.strptime(config['dates']['start'], '%Y-%m-%d %H:%M:%S') 
    end_date = datetime.strptime(config['dates']['end'], '%Y-%m-%d %H:%M:%S') 
    
    base_path = "/discover/nobackup/projects/caturne4/"
    experiment_name = "geos_cam_eval"

    aeronet_path = config['paths']['aeronet']
    forecast_path = config['paths']['forecast']
    FPanalysis_path = config['paths']['fp_forecast']

    
    output_dir = os.path.join(base_path, experiment_name, "comparison")
    min_points = 10   # of required hourly data points to be included in high-level stats
    ts_freq = config['temporal_resolution']
    evaluation = config['evaluation']
   

    
    main(
        start_date=start_date,
        end_date=end_date,
        base_path=base_path,
        experiment_name=experiment_name,
        output_dir=output_dir, 
        min_points=min_points, 
        ts_freq=ts_freq
    )
