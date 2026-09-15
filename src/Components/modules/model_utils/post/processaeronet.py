#!/usr/bin/env python3

'''
This code subsamples multiple GEOS model experiments according to hourly means of AERONET observations for both solar and lunar retrievals.
The output is a CSV file with high-level statistics for all sites with observations for the 
requested date range including means, standard deviations, RMSE, and correlations as well as
a time series CSV file at the requested temporal frequency (hourly, daily, monthly).
Since this utilizes aeronet.py from GMAOpyobs, the python path must be set first 
(example: export PYTHONPATH="/discover/nobackup/acollow/GMAOpyobs/src/:$PYTHONPATH")
'''

import os
import glob
import numpy as np
import pandas as pd
import xarray as xr
from datetime import datetime, timedelta
import warnings
import multiprocessing as mp
from functools import partial
import time
from pyobs.aeronet import AERONET_L2, VARS 

warnings.filterwarnings('ignore')

def parse_aeronet_file(file_info, start_date, end_date):
    """Reads an AERONET file (solar or lunar) using AERONET_L2 class, applies QC, and returns hourly means per station."""
    filepath, obs_type = file_info
    try:
        custom_vars = list(VARS) + ['440_870_Angstrom_Exponent']
        
        # Using standard AERONET V3 parser for both solar and lunar lev15 data
        aero = AERONET_L2(filepath, version=3, Vars=custom_vars, Verbose=False)
        
        if aero.nobs == 0:
            return None, "No data loaded"

        df = pd.DataFrame({
            'DateTime': aero.tyme,
            'lat': aero.Latitude,
            'lon': aero.Longitude,
            'station': aero.Location, 
            'aod_550': aero.AOT_550,
            'angstrom': getattr(aero, '440_870_Angstrom_Exponent', np.nan),
            'obs_type': obs_type
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
        
        hourly = df.groupby(['station', 'hour', 'obs_type']).agg({
            'aod_550': 'mean', 
            'angstrom': 'mean',
            'lat': 'first',
            'lon': 'first'
        }).reset_index()
        
        return hourly, "Success"
        
    except ValueError as e:
        return None, f"Parsing Error: {str(e)}"
    except Exception as e:
        return None, f"Exception: {type(e).__name__} - {str(e)}"


def processmodel(date_tuple, model_path_template):
    """Processes all stations for a single day against the GEOS model data"""
    date, daily_df = date_tuple
    
    for hour_val, group in daily_df.groupby('hour'):
        YYYY = hour_val.strftime('%Y')
        MM = hour_val.strftime('%m')
        YYYYMM = hour_val.strftime('%Y%m')
        YYYYMMDD = hour_val.strftime('%Y%m%d')
        HH = hour_val.strftime('%H')
        
        model_file = model_path_template.format(YYYY=YYYY, MM=MM, YYYYMM=YYYYMM, YYYYMMDD=YYYYMMDD, HH=HH)
        
        if os.path.exists(model_file):
            with xr.open_dataset(model_file) as ds_mod:
                for idx, row in group.iterrows():
                    mlon = row['lon'] + 360 if (ds_mod.lon.min() >= 0 and row['lon'] < 0) else row['lon']
                    
                    pt = ds_mod.sel(lat=row['lat'], lon=mlon, method='nearest')
                    
                    if 'time' in pt.dims:
                        pt = pt.isel(time=0)
                        
                    daily_df.at[idx, 'model_aod'] = float(pt['TOTEXTTAU'].values)
                    daily_df.at[idx, 'model_ang'] = float(pt['TOTANGSTR'].values)

    return daily_df


def main(start_date, end_date, model_base_path, experiment_names, min_points=50, ts_freq='none', output_base_path='./'):
    start_time = time.time()
    n_procs = max(1, mp.cpu_count() - 1)
    
    aeronet_solar_dir = "/css/gmao/dp/gds/AeroObs/AERONET.v3/Level1_5"
    aeronet_lunar_dir = "/discover/nobackup/acollow/aeroeval/aeronet_lunar/AOD_LUNAR/AOD15/ALL_POINTS"
    
    print(f" Locating AERONET solar & lunar files between {start_date.strftime('%Y-%m-%d')} and {end_date.strftime('%Y-%m-%d')}...")
    
    files_to_process = []
    
    # 1. Locate Solar Files
    current_date = start_date.replace(hour=0, minute=0, second=0)
    while current_date <= end_date:
        YYYY = current_date.strftime('%Y')
        MM = current_date.strftime('%m')
        YYYYMMDD = current_date.strftime('%Y%m%d')
        
        possible_filenames = [
            f"aeronet_v30.{YYYYMMDD}.txt",
            f"aeronet_v3.aod.{YYYYMMDD}.txt"
        ]
        
        for fname in possible_filenames:
            file_path = os.path.join(aeronet_solar_dir, f"Y{YYYY}", f"M{MM}", fname)
            
            if os.path.exists(file_path):
                files_to_process.append((file_path, 'solar'))
                break 
        
        current_date += timedelta(days=1)
        
    # 2. Locate Lunar Files
    lunar_files = glob.glob(os.path.join(aeronet_lunar_dir, "20130101_20260912_*.lev15"))
    for lf in lunar_files:
        files_to_process.append((lf, 'lunar'))
    
    if not files_to_process:
        print("No AERONET files found to process!")
        return
        
    print(f"Found {len(files_to_process)} total AERONET files to check. Parsing...")
        
    parse_func = partial(parse_aeronet_file, start_date=start_date, end_date=end_date)
    dfs = []
    error_counts = {}
    
    with mp.Pool(n_procs) as pool:
        for df, msg in pool.imap_unordered(parse_func, files_to_process):
            if df is not None:
                dfs.append(df)
            else:
                error_counts[msg] = error_counts.get(msg, 0) + 1

    print("\n--- AERONET Parse Summary ---")
    print(f"Successfully loaded data for {len(dfs)} file chunks.")
    
    if error_counts:
        print("\nReasons for skipped files (mostly expected due to out of date range for station files):")
        # Print top 5 reasons to not clog the terminal
        for reason, count in sorted(error_counts.items(), key=lambda x: x[1], reverse=True)[:5]: 
            print(f"  - {reason}: {count} files")
            
    if not dfs:
        print("\nNo valid data parsed. Exiting.")
        return
        
    master_df = pd.concat(dfs, ignore_index=True)
    master_df['date'] = master_df['hour'].dt.floor('D')
    
    print(f"\nTotal valid AERONET hours across all stations:")
    print(master_df['obs_type'].value_counts())
    
    for experiment_name in experiment_names:
        print(f"\n========================================================")
        print(f"Sampling Experiment: {experiment_name}")
        print(f"========================================================")
        
        output_dir = os.path.join(output_base_path, experiment_name, "AERONETsampled")
        os.makedirs(output_dir, exist_ok=True)
        
        # New model path logic: /discover/nobackup/projects/gmao/geos_fp_arch/f5430_fp/chem/Y2026/M06/f5430_fp.inst1_2d_hwl_Nx.20260601_0000z.nc4
        model_path_template = os.path.join(
            model_base_path, 
            experiment_name, 
            "chem", 
            "Y{YYYY}", 
            "M{MM}", 
            f"{experiment_name}.inst1_2d_hwl_Nx.{{YYYYMMDD}}_{{HH}}00z.nc4"
        )
        
        exp_df = master_df.copy()
        exp_df['model_aod'] = np.nan
        exp_df['model_ang'] = np.nan
        
        grouped = [(date, group) for date, group in exp_df.groupby('date')]
        process_func = partial(processmodel, model_path_template=model_path_template)
        processed_dfs = []
        
        with mp.Pool(n_procs) as pool:
            for i, result_df in enumerate(pool.imap_unordered(process_func, grouped)):
                processed_dfs.append(result_df)
                if (i+1) % 10 == 0 or (i+1) == len(grouped):
                    print(f"Processed {i+1}/{len(grouped)} days for {experiment_name}...")
                    
        final_df = pd.concat(processed_dfs, ignore_index=True)
        final_df = final_df.dropna(subset=['model_aod', 'model_ang'])
        
        # Filter by min_points ensuring station and obs_type grouping
        station_counts = final_df.groupby(['station', 'obs_type']).size()
        valid_groups = station_counts[station_counts >= min_points].reset_index()
        final_df = final_df.merge(valid_groups[['station', 'obs_type']], on=['station', 'obs_type'], how='inner')
        
        print(f"\nCalculating Statistics... ({len(final_df)} matched points)")
        stats = []
        for (station, obs_type), grp in final_df.groupby(['station', 'obs_type']):
            stats.append({
                'Experiment': experiment_name,
                'Station': station,
                'Obs_Type': obs_type,
                'Lat': grp['lat'].iloc[0],
                'Lon': grp['lon'].iloc[0],
                'N_Points': len(grp),
                
                'AERONET_AOD_Mean': grp['aod_550'].mean(),
                'AERONET_AOD_Std': grp['aod_550'].std(),
                'MODEL_AOD_Mean': grp['model_aod'].mean(),
                'MODEL_AOD_Std': grp['model_aod'].std(),
                
                'AERONET_Angstrom_Mean': grp['angstrom'].mean(),
                'AERONET_Angstrom_Std': grp['angstrom'].std(),
                'MODEL_Angstrom_Mean': grp['model_ang'].mean(),
                'MODEL_Angstrom_Std': grp['model_ang'].std(),
                
                'MODEL_AOD_RMSE': np.sqrt(((grp['model_aod'] - grp['aod_550'])**2).mean()),
                'MODEL_Angstrom_RMSE': np.sqrt(((grp['model_ang'] - grp['angstrom'])**2).mean()),
                
                'MODEL_AOD_Correlation': grp['model_aod'].corr(grp['aod_550']),
                'MODEL_Angstrom_Correlation': grp['model_ang'].corr(grp['angstrom'])
            })
            
        date_str = f"{start_date.strftime('%Y%m%d')}_{end_date.strftime('%Y%m%d')}"
        
        if stats:
            out_df = pd.DataFrame(stats)  
            output_file = os.path.join(output_dir, f"{experiment_name}_aeronet_comparison_stats_{date_str}.csv")
            out_df.to_csv(output_file, index=False)
            print(f"Stats saved to {output_file}")
        else:
            print(f"No stations met the minimum threshold of {min_points} valid matched points to calculate statistics.")

        if ts_freq != 'none' and not final_df.empty:
            print(f"Generating {ts_freq} timeseries data...")
            
            cols_to_keep = ['hour', 'station', 'obs_type', 'lat', 'lon', 'aod_550', 'model_aod', 'angstrom', 'model_ang']
            ts_df = final_df[cols_to_keep].copy()
            ts_df.rename(columns={'hour': 'time', 'aod_550': 'aeronet_aod', 'angstrom': 'aeronet_angstrom'}, inplace=True)
            
            ts_df.insert(0, 'experiment', experiment_name)
            
            if ts_freq == 'daily':
                ts_df['time'] = ts_df['time'].dt.floor('D')
                ts_df = ts_df.groupby(['experiment', 'station', 'obs_type', 'time', 'lat', 'lon']).mean(numeric_only=True).reset_index()
            elif ts_freq == 'monthly':
                ts_df['time'] = ts_df['time'].dt.to_period('M').dt.to_timestamp()
                ts_df = ts_df.groupby(['experiment', 'station', 'obs_type', 'time', 'lat', 'lon']).mean(numeric_only=True).reset_index()
                
            ts_output_file = os.path.join(output_dir, f"{experiment_name}_aeronet_timeseries_{ts_freq}_{date_str}.csv")
            
            ts_df.to_csv(ts_output_file, index=False)
            print(f"Time series data saved to {ts_output_file}")

    print(f"\nTotal script time: {(time.time() - start_time)/60:.1f} minutes!")


if __name__ == "__main__":
    
    start_date = datetime(2026, 6, 1, 0, 0, 0)
    end_date = datetime(2026, 8, 31, 23, 59, 59)
    
    # Updated Base Path to GEOS_FP
    model_base_path = "/discover/nobackup/projects/gmao/geos_fp_arch"
    
    # Path where outputs will be stored locally to avoid archive permission errors
    output_base_path = os.getcwd()
    
    experiment_names = [
        "f5430_fp"
    ]
    
    min_points = 4   
    ts_freq = 'daily' 
    
    main(
        start_date=start_date,
        end_date=end_date,
        model_base_path=model_base_path,
        experiment_names=experiment_names,
        min_points=min_points, 
        ts_freq=ts_freq,
        output_base_path=output_base_path
    )
