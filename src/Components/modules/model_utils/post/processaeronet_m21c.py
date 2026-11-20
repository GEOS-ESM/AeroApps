'''
This code subsamples MERRA-2 and MERRA-21C according to hourly means of AERONET observations.
The output is a CSV file with high-level statistics for all sites with observations for the 
requested year including means, standard deviations, RMSE, and correlations.
It also outputs a time series CSV file at the requested temporal frequency using the --ts-freq flag.

To Do:
-consider expanding to additional wavelengths
-remove hardwired paths for MERRA-2* and replace with flexible experiment paths
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

def parse_aeronet_file(filepath, years_to_process):
    """Reads a single AERONET file, computes AOD 550 using 440-675 AE in needed, and returns hourly mean AOD 550 and 440-870 AE for model comparison."""
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
            return None, f"Missing Date/Time columns in {station}"
            
        df['DateTime'] = pd.to_datetime(df['Date(dd:mm:yyyy)'] + ' ' + df['Time(hh:mm:ss)'], format='%d:%m:%Y %H:%M:%S', errors='coerce')
        df = df.dropna(subset=['DateTime'])
        
        df = df[df['DateTime'].dt.year.isin(years_to_process)]
        if df.empty: 
            return None, "No data in specified years"
        
        lat = float(df['Site_Latitude(Degrees)'].iloc[0])
        lon = float(df['Site_Longitude(Degrees)'].iloc[0])
        
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

# -------------------------------------------------------------------
# STAGE 2: Process Models by DAY
# -------------------------------------------------------------------
def process_day(date_tuple, merra2_dir, merra21c_dir):
    """Processes all stations for a single day, opening files only once."""
    date, daily_df = date_tuple
    y, m, d = date.year, date.month, date.day
    
    m2_file = os.path.join(merra2_dir, f"Y{y}", f"M{m:02d}", f"MERRA2.tavg1_2d_aer_Nx.{y}{m:02d}{d:02d}.nc4")
    if os.path.exists(m2_file):
        try:
            with xr.open_dataset(m2_file) as ds_m2:
                for idx, row in daily_df.iterrows():
                    mlon = row['lon'] + 360 if (ds_m2.lon.min() >= 0 and row['lon'] < 0) else row['lon']
                    try:
                        pt = ds_m2.sel(time=row['hour'], lat=row['lat'], lon=mlon, method='nearest')
                        daily_df.at[idx, 'm2_aod'] = float(pt['TOTEXTTAU'].values)
                        daily_df.at[idx, 'm2_ang'] = float(pt['TOTANGSTR'].values)
                    except Exception:
                        pass 
        except Exception:
            pass 

    for hour_val, group in daily_df.groupby('hour'):
        h = hour_val.hour
        m21c_file = os.path.join(merra21c_dir, f"Y{y}", f"M{m:02d}", f"e5303_m21c_jan18.aer_inst_1hr_glo_L1152x721_slv.{y}-{m:02d}-{d:02d}T{h:02d}00Z.nc4")
        
        if os.path.exists(m21c_file):
            try:
                with xr.open_dataset(m21c_file) as ds_m21c:
                    for idx, row in group.iterrows():
                        mlon = row['lon'] + 360 if (ds_m21c.lon.min() >= 0 and row['lon'] < 0) else row['lon']
                        try:
                            pt = ds_m21c.sel(lat=row['lat'], lon=mlon, method='nearest')
                            daily_df.at[idx, 'm21c_aod'] = float(pt['TOTEXTTAU'].values)
                            daily_df.at[idx, 'm21c_ang'] = float(pt['TOTANGSTR'].values)
                        except Exception:
                            pass
            except Exception:
                pass 

    return daily_df

# -------------------------------------------------------------------
# MAIN EXECUTION
# -------------------------------------------------------------------
def main(years_to_process=[2018, 2019, 2020, 2021, 2022, 2023, 2024], output_dir="./aeronet_merra_comparison/", min_points=50, ts_freq='none'):
    start_time = time.time()
    n_procs = max(1, mp.cpu_count() - 1)
    os.makedirs(output_dir, exist_ok=True)
    
    aeronet_dir = "/discover/nobackup/acollow/aeroeval/aeronet/m21c/AOD/AOD20/ALL_POINTS/"
    merra2_dir = "/discover/nobackup/projects/gmao/merra2/data/pub/products/MERRA2_all/"
    merra21c_dir = "/discover/nobackup/projects/gmao/merra21c/archive/e5303_m21c_jan18/chem/"
    
    print("STAGE 1: Parsing AERONET files...")
    files = glob.glob(os.path.join(aeronet_dir, "19930101_20260411_*.lev20"))
    
    if not files:
        print(f"CRITICAL ERROR: No files found in {aeronet_dir}")
        return
        
    parse_func = partial(parse_aeronet_file, years_to_process=years_to_process)
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
    
    if not dfs:
        print("\nCRITICAL ERROR: No dataframes to process. Exiting.")
        return
        
    master_df = pd.concat(dfs, ignore_index=True)
    master_df['date'] = master_df['hour'].dt.floor('D')
    print(f"\nTotal valid AERONET hours across all stations: {len(master_df)}")
    
    print("\nSTAGE 2: Extracting Model Data (Time-based)...")
    master_df['m2_aod'] = np.nan
    master_df['m2_ang'] = np.nan
    master_df['m21c_aod'] = np.nan
    master_df['m21c_ang'] = np.nan
    
    grouped = [(date, group) for date, group in master_df.groupby('date')]
    process_func = partial(process_day, merra2_dir=merra2_dir, merra21c_dir=merra21c_dir)
    processed_dfs = []
    
    with mp.Pool(n_procs) as pool:
        for i, result_df in enumerate(pool.imap_unordered(process_func, grouped)):
            processed_dfs.append(result_df)
            if (i+1) % 50 == 0 or (i+1) == len(grouped):
                print(f"Processed {i+1}/{len(grouped)} days...")
                
    final_df = pd.concat(processed_dfs, ignore_index=True)
    
    # Clean model data using safe ranges
    final_df.loc[(final_df['m2_aod'] < 0) | (final_df['m2_aod'] >= 10), 'm2_aod'] = np.nan
    final_df.loc[(final_df['m21c_aod'] < 0) | (final_df['m21c_aod'] >= 10), 'm21c_aod'] = np.nan
    final_df.loc[(final_df['m2_ang'] < -1) | (final_df['m2_ang'] > 3), 'm2_ang'] = np.nan
    final_df.loc[(final_df['m21c_ang'] < -1) | (final_df['m21c_ang'] > 3), 'm21c_ang'] = np.nan
    
    # Keep only instances where we have data for AERONET, MERRA-2, and MERRA-21C simultaneously
    final_df = final_df.dropna(subset=['m2_aod', 'm21c_aod', 'm2_ang', 'm21c_ang'])
    
    # Filter for stations meeting the min_points threshold
    valid_stations = final_df['station'].value_counts()[final_df['station'].value_counts() >= min_points].index
    final_df = final_df[final_df['station'].isin(valid_stations)].copy()
    
    print(f"\nSTAGE 3: Calculating Statistics... ({len(final_df)} matched hour-station points from {len(valid_stations)} stations)")
    stats = []
    for station, grp in final_df.groupby('station'):
        stats.append({
            'Station': station,
            'Lat': grp['lat'].iloc[0],
            'Lon': grp['lon'].iloc[0],
            'N_Points': len(grp),
            
            'AERONET_AOD_Mean': grp['aod_550'].mean(),
            'AERONET_AOD_Std': grp['aod_550'].std(),
            'MERRA2_AOD_Mean': grp['m2_aod'].mean(),
            'MERRA2_AOD_Std': grp['m2_aod'].std(),
            'MERRA21C_AOD_Mean': grp['m21c_aod'].mean(),
            'MERRA21C_AOD_Std': grp['m21c_aod'].std(),
            
            'AERONET_Angstrom_Mean': grp['angstrom'].mean(),
            'AERONET_Angstrom_Std': grp['angstrom'].std(),
            'MERRA2_Angstrom_Mean': grp['m2_ang'].mean(),
            'MERRA2_Angstrom_Std': grp['m2_ang'].std(),
            'MERRA21C_Angstrom_Mean': grp['m21c_ang'].mean(),
            'MERRA21C_Angstrom_Std': grp['m21c_ang'].std(),
            
            'MERRA2_AOD_RMSE': np.sqrt(((grp['m2_aod'] - grp['aod_550'])**2).mean()),
            'MERRA21C_AOD_RMSE': np.sqrt(((grp['m21c_aod'] - grp['aod_550'])**2).mean()),
            'MERRA2_Angstrom_RMSE': np.sqrt(((grp['m2_ang'] - grp['angstrom'])**2).mean()),
            'MERRA21C_Angstrom_RMSE': np.sqrt(((grp['m21c_ang'] - grp['angstrom'])**2).mean()),
            
            'MERRA2_AOD_Correlation': grp['m2_aod'].corr(grp['aod_550']),
            'MERRA21C_AOD_Correlation': grp['m21c_aod'].corr(grp['aod_550']),
            'MERRA2_Angstrom_Correlation': grp['m2_ang'].corr(grp['angstrom']),
            'MERRA21C_Angstrom_Correlation': grp['m21c_ang'].corr(grp['angstrom'])
        })
        
    year_str = f"{years_to_process[0]}_{years_to_process[-1]}" if len(years_to_process) > 1 else str(years_to_process[0])
    
    if stats:
        out_df = pd.DataFrame(stats).round(4)
        output_file = os.path.join(output_dir, f"aeronet_merra_comparison_stats_{year_str}.csv")
        out_df.to_csv(output_file, index=False)
        print(f"Stats saved to {output_file}")
    else:
        print(f"No stations met the minimum threshold of {min_points} valid matched points to calculate statistics.")

    # -------------------------------------------------------------------
    # STAGE 4: Output Time Series Data
    # -------------------------------------------------------------------
    if ts_freq != 'none' and not final_df.empty:
        print(f"\nSTAGE 4: Generating {ts_freq} timeseries data...")
        
        # Select and rename columns for clarity
        cols_to_keep = ['hour', 'station', 'lat', 'lon', 'aod_550', 'm2_aod', 'm21c_aod', 'angstrom', 'm2_ang', 'm21c_ang']
        ts_df = final_df[cols_to_keep].copy()
        ts_df.rename(columns={'hour': 'time', 'aod_550': 'aeronet_aod', 'angstrom': 'aeronet_angstrom'}, inplace=True)
        
        # Apply Temporal Aggregation
        if ts_freq == 'daily':
            ts_df['time'] = ts_df['time'].dt.floor('D')
            ts_df = ts_df.groupby(['station', 'time', 'lat', 'lon']).mean().reset_index()
        elif ts_freq == 'monthly':
            ts_df['time'] = ts_df['time'].dt.to_period('M').dt.to_timestamp()
            ts_df = ts_df.groupby(['station', 'time', 'lat', 'lon']).mean().reset_index()
            
        ts_output_file = os.path.join(output_dir, f"aeronet_merra_timeseries_{ts_freq}_{year_str}.csv")
        
        # Save roundings for cleaner files
        cols_to_round = ['aeronet_aod', 'm2_aod', 'm21c_aod', 'aeronet_angstrom', 'm2_ang', 'm21c_ang']
        ts_df[cols_to_round] = ts_df[cols_to_round].round(4)
        
        ts_df.to_csv(ts_output_file, index=False)
        print(f"Time series data saved to {ts_output_file}")

    print(f"Total time: {(time.time() - start_time)/60:.1f} minutes!")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--years', nargs='+', type=int, default=[2018, 2019, 2020, 2021, 2022, 2023, 2024])
    parser.add_argument('--output-dir', type=str, default="./aeronet_merra_comparison/")
    parser.add_argument('--min-points', type=int, default=50, help="Minimum valid matched hourly points required to include a station in the final statistics.")
    parser.add_argument('--ts-freq', type=str, choices=['none', 'hourly', 'daily', 'monthly'], default='none', help="Temporal frequency for time series output file.")
    args = parser.parse_args()
    
    main(years_to_process=args.years, output_dir=args.output_dir, min_points=args.min_points, ts_freq=args.ts_freq)
