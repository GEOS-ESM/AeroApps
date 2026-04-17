#!/usr/bin/env python3

'''
This code subsamples a single GEOS model experiment according to hourly means of AERONET observations.
The output is a CSV file with high-level statistics for all sites with observations for the 
requested date range including means, standard deviations, RMSE, and correlations as well as
a time series CSV file at the requested temporal frequency (hourly, daily, monthly).
'''

import os
import glob
import numpy as np
import pandas as pd
import xarray as xr
from datetime import datetime
import warnings
import multiprocessing as mp
from functools import partial
import time
warnings.filterwarnings('ignore')

def parse_aeronet_file(filepath, start_date, end_date):
    """Reads a single AERONET file, computes AOD 550, applies QC, and returns hourly mean AOD 550 and 440-870 AE."""
    station = os.path.basename(filepath).split("_", 2)[2].replace(".lev20", "")
    
    try:
        try:
            with open(filepath, 'r', encoding='utf-8') as f:
                lines = f.readlines()
                header_line = lines[6].strip()
            df = pd.read_csv(filepath, skiprows=7, header=None, sep=',', encoding='utf-8')
        except UnicodeDecodeError:
            with open(filepath, 'r', encoding='latin-1') as f:
                lines = f.readlines()
                header_line = lines[6].strip()
            df = pd.read_csv(filepath, skiprows=7, header=None, sep=',', encoding='latin-1')
            
        column_names = header_line.split(',')
        
        if len(column_names) > len(df.columns):
            column_names = column_names[:len(df.columns)]
        elif len(df.columns) > len(column_names):
            column_names.extend([f"Unnamed_{i}" for i in range(len(df.columns) - len(column_names))])
            
        df.columns = column_names
        df = df.loc[:, ~df.columns.duplicated()].copy()
            
        if 'Date(dd:mm:yyyy)' not in df.columns or 'Time(hh:mm:ss)' not in df.columns:
            return None, f"Missing Date/Time columns"
            
        df['DateTime'] = pd.to_datetime(df['Date(dd:mm:yyyy)'] + ' ' + df['Time(hh:mm:ss)'], format='%d:%m:%Y %H:%M:%S', errors='coerce')
        df = df.dropna(subset=['DateTime'])
        
        # Filter by specific date range
        df = df[(df['DateTime'] >= start_date) & (df['DateTime'] <= end_date)]
        if df.empty: 
            return None, "No data in specified date range"
        
        lat = float(df['Site_Latitude(Degrees)'].iloc[0])
        lon = float(df['Site_Longitude(Degrees)'].iloc[0])
        
        # QUALITY CONTROL: Angstrom Exponents between -1 and 3
        ang_440_675 = None
        if '440-675_Angstrom_Exponent' in df.columns:
            df['440-675_Angstrom_Exponent'] = pd.to_numeric(df['440-675_Angstrom_Exponent'], errors='coerce')
            valid_mask_675 = (df['440-675_Angstrom_Exponent'] >= -1) & (df['440-675_Angstrom_Exponent'] <= 3)
            ang_440_675 = df['440-675_Angstrom_Exponent'].where(valid_mask_675)
            
        ang_440_870 = None
        if '440-870_Angstrom_Exponent' in df.columns:
            df['440-870_Angstrom_Exponent'] = pd.to_numeric(df['440-870_Angstrom_Exponent'], errors='coerce')
            valid_mask_870 = (df['440-870_Angstrom_Exponent'] >= -1) & (df['440-870_Angstrom_Exponent'] <= 3)
            ang_440_870 = df['440-870_Angstrom_Exponent'].where(valid_mask_870)

        if ang_440_870 is None or ang_440_870.isna().all():
            return None, "Missing valid 440-870 Angstrom Exponent"

        # QUALITY CONTROL: AOD between 0 and 10
        aod_columns = ['AOD_551nm', 'AOD_550nm', 'AOD_532nm', 'AOD_531nm', 'AOD_555nm', 'AOD_560nm']
        aod_550 = None
        
        for col in aod_columns:
            if col in df.columns:
                df[col] = pd.to_numeric(df[col], errors='coerce')
                valid_aod = df[col][(df[col] >= 0) & (df[col] < 10)]
                if len(valid_aod) > 0:
                    aod_550 = df[col].copy()
                    aod_550[(df[col] < 0) | (df[col] >= 10)] = np.nan
                    break
        
        # Fallback to computing AOD at 550nm using 440nm and 440-675 AE
        if aod_550 is None:
            if 'AOD_440nm' in df.columns and ang_440_675 is not None:
                df['AOD_440nm'] = pd.to_numeric(df['AOD_440nm'], errors='coerce')
                aod_440 = df['AOD_440nm']
                
                valid_mask = (aod_440 >= 0) & (aod_440 < 10) & (~aod_440.isna()) & (~ang_440_675.isna())
                
                if valid_mask.sum() > 0:
                    aod_550 = pd.Series(np.nan, index=df.index)
                    aod_550[valid_mask] = aod_440[valid_mask] * (0.550 / 0.440) ** (-ang_440_675[valid_mask])
            
        if aod_550 is None or aod_550.isna().all():
            return None, "Cannot obtain or compute AOD at 550nm"
            
        df['aod_550'] = aod_550
        df['angstrom'] = ang_440_870
        df['hour'] = df['DateTime'].dt.floor('H')
        
        hourly = df.groupby('hour').agg({'aod_550':'mean', 'angstrom':'mean'}).dropna().reset_index()
        
        if hourly.empty:
            return None, "No valid paired AOD/Angstrom data"
            
        hourly['station'] = station
        hourly['lat'] = lat
        hourly['lon'] = lon
        
        return hourly, "Success"
        
    except Exception as e:
        return None, f"Exception: {type(e).__name__} - {str(e)}"

def processmodel(date_tuple, model_path_template):
    """Processes all stations for a single day"""
    date, daily_df = date_tuple
    
    for hour_val, group in daily_df.groupby('hour'):
        YYYYMM = hour_val.strftime('%Y%m')
        YYYYMMDD = hour_val.strftime('%Y%m%d')
        HH = hour_val.strftime('%H')
        model_file = model_path_template.format(YYYYMM=YYYYMM, YYYYMMDD=YYYYMMDD, HH=HH)
        
        if os.path.exists(model_file):
            with xr.open_dataset(model_file) as ds_mod:
                for idx, row in group.iterrows():
                    # Handle longitude wrapping if model is 0-360 and AERONET is -180 to 180
                    mlon = row['lon'] + 360 if (ds_mod.lon.min() >= 0 and row['lon'] < 0) else row['lon']
                    
                    pt = ds_mod.sel(lat=row['lat'], lon=mlon, method='nearest')
                    
                    if 'time' in pt.dims:
                        pt = pt.isel(time=0)
                        
                    daily_df.at[idx, 'model_aod'] = float(pt['TOTEXTTAU'].values)
                    daily_df.at[idx, 'model_ang'] = float(pt['TOTANGSTR'].values)

    return daily_df

def main(start_date, end_date, base_path, experiment_name, output_dir="./aeronet_model_comparison/", min_points=50, ts_freq='none'):
    start_time = time.time()
    n_procs = max(1, mp.cpu_count() - 1)
    os.makedirs(output_dir, exist_ok=True)
    
    # Generate the model path template dynamically
    model_path_template = os.path.join(base_path, experiment_name, "holding", "inst2d_hwl_x", "{YYYYMM}", f"{experiment_name}.inst2d_hwl_x.{{YYYYMMDD}}_{{HH}}00z.nc4")
    
    aeronet_dir = "/discover/nobackup/acollow/aeroeval/aeronet/m21c/AOD/AOD20/ALL_POINTS/"
    
    print(f" Parsing AERONET files between {start_date.strftime('%Y-%m-%d')} and {end_date.strftime('%Y-%m-%d')}...")
    files = glob.glob(os.path.join(aeronet_dir, "*.lev20"))
    
    if not files:
        print(f"No files found in {aeronet_dir}")
        return
        
    parse_func = partial(parse_aeronet_file, start_date=start_date, end_date=end_date)
    dfs = []
    error_counts = {}
    
    with mp.Pool(n_procs) as pool:
        for df, msg in pool.imap_unordered(parse_func, files):
            if df is not None:
                dfs.append(df)
            else:
                error_counts[msg] = error_counts.get(msg, 0) + 1

    print("\n--- AERONET Parse Summary ---")
    print(f"Successfully loaded {len(dfs)} stations.")
    
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
    
    grouped = [(date, group) for date, group in master_df.groupby('date')]
    process_func = partial(processmodel, model_path_template=model_path_template)
    processed_dfs = []
    
    with mp.Pool(n_procs) as pool:
        for i, result_df in enumerate(pool.imap_unordered(process_func, grouped)):
            processed_dfs.append(result_df)
            if (i+1) % 10 == 0 or (i+1) == len(grouped):
                print(f"Processed {i+1}/{len(grouped)} days...")
                
    final_df = pd.concat(processed_dfs, ignore_index=True)
    
    # Sample for available AERONET Obs
    final_df = final_df.dropna(subset=['model_aod', 'model_ang'])
    
    # Filter for stations meeting the min_points threshold
    valid_stations = final_df['station'].value_counts()[final_df['station'].value_counts() >= min_points].index
    final_df = final_df[final_df['station'].isin(valid_stations)].copy()
    
    print(f"\nCalculating Statistics... ({len(final_df)} matched hour-station points from {len(valid_stations)} stations)")
    stats = []
    for station, grp in final_df.groupby('station'):
        stats.append({
            'Experiment': experiment_name,
            'Station': station,
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
        out_df = pd.DataFrame(stats)  # Rounding removed
        output_file = os.path.join(output_dir, f"{experiment_name}_aeronet_comparison_stats_{date_str}.csv")
        out_df.to_csv(output_file, index=False)
        print(f"Stats saved to {output_file}")
    else:
        print(f"No stations met the minimum threshold of {min_points} valid matched points to calculate statistics.")

    if ts_freq != 'none' and not final_df.empty:
        print(f"\nGenerating {ts_freq} timeseries data...")
        
        # Select and rename columns for clarity
        cols_to_keep = ['hour', 'station', 'lat', 'lon', 'aod_550', 'model_aod', 'angstrom', 'model_ang']
        ts_df = final_df[cols_to_keep].copy()
        ts_df.rename(columns={'hour': 'time', 'aod_550': 'aeronet_aod', 'angstrom': 'aeronet_angstrom'}, inplace=True)
        
        ts_df.insert(0, 'experiment', experiment_name)
        
        # Apply Temporal Aggregation
        if ts_freq == 'daily':
            ts_df['time'] = ts_df['time'].dt.floor('D')
            ts_df = ts_df.groupby(['experiment', 'station', 'time', 'lat', 'lon']).mean().reset_index()
        elif ts_freq == 'monthly':
            ts_df['time'] = ts_df['time'].dt.to_period('M').dt.to_timestamp()
            ts_df = ts_df.groupby(['experiment', 'station', 'time', 'lat', 'lon']).mean().reset_index()
            
        ts_output_file = os.path.join(output_dir, f"{experiment_name}_aeronet_timeseries_{ts_freq}_{date_str}.csv")
        
        ts_df.to_csv(ts_output_file, index=False)
        print(f"Time series data saved to {ts_output_file}")

    print(f"Total time: {(time.time() - start_time)/60:.1f} minutes!")

if __name__ == "__main__":
    
    start_date = datetime(2024, 1, 1, 0, 0, 0)
    end_date = datetime(2024, 1, 31, 23, 59, 59)
    
    base_path = "/discover/nobackup/projects/gmao/geos_aerosols/acollow/"
    experiment_name = "c180R_qfed3-2_scaled"
    
    # Optional Output configuration
    output_dir = os.path.join(base_path, experiment_name, "AERONETsampled")
    min_points = 10  # of required hourly data points to be included in high-level stats
    ts_freq = 'daily' # Set to 'hourly', 'daily', or 'monthly' 
    
    # =========================================================================
    
    main(
        start_date=start_date,
        end_date=end_date,
        base_path=base_path,
        experiment_name=experiment_name,
        output_dir=output_dir, 
        min_points=min_points, 
        ts_freq=ts_freq
    )
