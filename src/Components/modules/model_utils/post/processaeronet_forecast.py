#!/usr/bin/env python3

'''
This code can sample two GEOS model datasets and calculates statistics against hourly means of AERONET observations.
Originally intended for GEOS-CAM and GEOS-FP so changes may be necessary for other datasets.
To cahnge settings including start and end dates, data directories, output directories, etc, see config.yaml.
The output consists of 4 .csv files and one netCDF file.
The .csv files consist of high level statistics for all sites with observations for the 
requested date range including means, standard deviations, RMSE, and correlations as well as
a time series CSV file at the requested temporal frequency (hourly, daily, monthly).
Since this utilizes aeronet.py from GMAOpyobs, the python path must be set first (export PYTHONPATH="/discover/nobackup/caturne4/AeroApps/Python")
'''
 
import os
import sys
sys.path.append('/discover/nobackup/caturne4/AeroApps/Python') # changes python path to allow for pyobs import
import numpy as np
import pandas as pd
import xarray as xr
from datetime import datetime, timedelta
import warnings
import multiprocessing as mp
from functools import partial
import time
from pyobs.aeronet import AERONET_L2, VARS 
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
        
        df = pd.DataFrame({               # creating dataframe from aeronet data
            'DateTime': aero.tyme,
            'lats': aero.Latitude,
            'lons': aero.Longitude,
            'station': aero.Location, 
            'aod_550': aero.AOT_550,
            'angstrom': getattr(aero, '440_870_Angstrom_Exponent', np.nan) 
        })

        # subsetting to start and end dates
        df = df[(df['DateTime'] >= start_date) & (df['DateTime'] <= end_date)] 
        if df.empty: 
            return None, "No data in specified date range"


        # cleaning the data
        df.replace(-999.0, np.nan, inplace=True)
        df.loc[(df['aod_550'] < 0) | (df['aod_550'] >= 10), 'aod_550'] = np.nan 
        valid_mask_870 = (df['angstrom'] >= -1) & (df['angstrom'] <= 3)
        df.loc[~valid_mask_870, 'angstrom'] = np.nan


        if df.empty:
            return None, "No valid paired AOD/Angstrom data after QC"

        # organizing the data
        df['obs hour'] = df['DateTime'].dt.floor('h')
        
        hourly = df.groupby(['station', 'obs hour']).agg({
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
    
    station_cache = {}

    for hour_val, group in daily_df.groupby('obs hour'):
        # creatting date variables
        YYYYMM = hour_val.strftime('%Y%m')
        YYYYMMDD = hour_val.strftime('%Y%m%d')
        MM = hour_val.strftime('%m')
        HH = hour_val.strftime('%H')
        HHmm = hour_val.strftime('%H%M')
       
        analysis_file = analysis_path_template.format(YYYYMM=YYYYMM, YYYYMMDD=YYYYMMDD, HH=HH, MM=MM, HHmm = HHmm)
        if os.path.exists(analysis_file):
            try:
                with xr.open_dataset(analysis_file, decode_times=False) as ds_mod:
                    for idx, row in group.iterrows():
                        station = row['station']
                        
                        if station not in station_cache:
                            target_lat = row['lats']
                            target_lon = row['lons']
                            try:
                                lon_diff = (ds_mod['lons'] - target_lon +180)%360-180
                                dist = (ds_mod['lats'] - target_lat)**2 + lon_diff**2
                                nearest_point = dist.compute().argmin(dim=['Xdim', 'Ydim'])
                                mlon = target_lon + 360 if (ds_mod.lons.min() >= 0 and target_lon < 0) else target_lon
                                station_cache[station] = {'point': nearest_point, 'mlon': mlon, 'method': 'argmin'}
                            except Exception:
                                mlon = target_lon + 360 if (ds_mod.lon.min() >= 0 and target_lon < 0) else target_lon
                                station_cache[station] = {'lat': target_lat, 'mlon': mlon, 'method': 'sel'}
                                
                        cache = station_cache[station]
                        
                        if cache['method'] == 'argmin':
                            pt = ds_mod.isel(cache['point'])
                        else:
                            pt = ds_mod.sel(lat=cache['lat'], lon=cache['mlon'], method='nearest')
    
                        try:
                            if 'time' in pt.dims:
                                pt = pt.isel(time=0)
                                
                            # gathering data for analysis
                            daily_df.at[idx, 'analysis_aod'] = float(pt['TOTEXTTAU'].values) 
                            daily_df.at[idx, 'analysis_ang'] = float(pt['TOTANGSTR'].values) 
                        except IndexError:
                            continue
            except PermissionError as e:
                print('Permissioin to file denied', e)
                continue
        else:
            print(f'analysis file {analysis_file} not found')

    return daily_df


#####################################################################################
# Process the model data
#####################################################################################
def processmodel(date_tuple, model_path_template):
    """Processes all stations for a single day against the GEOS model data"""
    date, daily_df = date_tuple
    processed_records = []
    
    station_cache = {}
    climate_cache = {}
    
    day_start = daily_df['obs hour'].min() - pd.Timedelta(hours=fcst_period)
    day_end = daily_df['obs hour'].max()
    day_start = max(day_start, start_date)
    day_end = min(day_end, end_date)       
    
    valid_init_times = pd.date_range(start=day_start, end=day_end, freq=time_step)

    for initiation_time in valid_init_times:
        # setting variables for model initalization time
        initiation_YYYYMMDD = datetime.strftime(initiation_time, '%Y%m%d')
        initiation_YYYY = datetime.strftime(initiation_time, '%Y')
        initiation_MM = datetime.strftime(initiation_time, '%m')
        initiation_DD = datetime.strftime(initiation_time, '%d')
        initiation_HH = datetime.strftime(initiation_time, '%H')
        
        for hour_val, group in daily_df.groupby('obs hour'):
            # setting variables for model forecast time
            YYYYMM = hour_val.strftime('%Y%m'), # model forecast year and month
            YYYYMMDD = hour_val.strftime('%Y%m%d') # model forecast year, month, and day
            MM = hour_val.strftime('%m') # model forecast month
            HHmm = hour_val.strftime('%H%M') # model forecast hour
            
            fcst_length = hour_val - initiation_time
            if fcst_length < pd.Timedelta(0):
                continue
            if fcst_length > pd.Timedelta(hours=fcst_period):
                continue
                
            try:
                model_file = model_path_template.format(initiation_YYYYMMDD = initiation_YYYYMMDD, initiation_HH = initiation_HH,
                                                        YYYYMM=YYYYMM, YYYYMMDD=YYYYMMDD, HH=initiation_HH, HHmm = HHmm)
            except KeyError:
                model_file = model_path_template.format(initiation_YYYYMMDD = initiation_YYYYMMDD,
                                                        initiation_YYYY = initiation_YYYY, 
                                                        initiation_MM = initiation_MM,
                                                        initiation_DD = initiation_DD,
                                                        initiation_HH = initiation_HH,
                                                        YYYYMM=YYYYMM, YYYYMMDD=YYYYMMDD, 
                                                        HH=initiation_HH, HHmm = HHmm)
                
        # Adding Climatology
            climate_path = '/discover/nobackup/acollow/aeroeval/geosit/monthly_climatologies/'
            climate_hour = (hour_val.hour // 6) * 6
            climate_file = f'{climate_path}geos_it_climatology_{climate_hour:02d}z_{MM}_2003_2022.nc'
            
            if os.path.exists(model_file) and os.path.exists(climate_file):
                if climate_file not in climate_cache:
                    climate_cache[climate_file] = xr.open_dataset(climate_file, decode_times = False)
                ds_climate = climate_cache[climate_file]
                with xr.open_dataset(model_file, decode_times=False) as ds_mod:
                    for idx, row in group.iterrows():
                        new_row = row.copy()
                        station = row['station']
                        
                        if station not in station_cache:
                            target_lat = row['lats']
                            target_lon = row['lons']
                            
                            try:
                                lon_diff = (ds_mod['lons'] - target_lon +180)%360-180
                                dist = (ds_mod['lats'] - target_lat)**2 + lon_diff**2
                                nearest_point = dist.compute().argmin(dim=['Xdim', 'Ydim'])
                                mlon = target_lon + 360 if (ds_mod.lons.min() >= 0 and target_lon < 0) else target_lon
                                
                                station_cache[station] = {'point': nearest_point, 'mlon': mlon, 'method': 'argmin'}
                            except Exception:
                                # print('in Exception\n~~~~~\n')
                                mlon = target_lon + 360 if (ds_mod.lon.min() >= 0 and target_lon < 0) else target_lon
                                station_cache[station] = {'lat': target_lat, 'mlon': mlon, 'method': 'sel'}

                        cache = station_cache[station]
                        
                        if cache['method'] == 'argmin':
                            pt = ds_mod.isel(cache['point'])
                            pt_climate = ds_climate.sel(lat=row['lats'], lon=cache['mlon'], method='nearest')
                        else:
                            pt = ds_mod.sel(lat=cache['lat'], lon=cache['mlon'], method='nearest')
                            pt_climate = ds_climate.sel(lat=cache['lat'], lon=cache['mlon'], method='nearest')

                        try: 
                            if 'time' in pt.dims:
                                pt = pt.isel(time=0)
    
                            # gathering values for stats
                            new_row['model_aod'] = float(pt['TOTEXTTAU'].squeeze().values.flatten()[0]) 
                            new_row['model_ang'] = float(pt['TOTANGSTR'].squeeze().values.flatten()[0])
                            new_row['climatology_aod'] = float(pt_climate['TOTEXTTAU'].squeeze().values.flatten()[0])
                            new_row['initialization time'] = initiation_time
                            new_row['initialization hour'] = initiation_HH
                            new_row['fcst length'] = fcst_length
                            
                            processed_records.append(new_row)
                        except IndexError:
                            continue

    
    
    if processed_records:
        return pd.DataFrame(processed_records)
    else:
        return pd.DataFrame(columns= list(daily_df.columns) + ['model_aod', 'model_ang', 'initialization time', 'initialization hour',  'fcst length'])
#####################################################################################
# Process the analysis data
#####################################################################################
def displayStats(data):
    # used to sort stats into seperate forecast hours

    stats = []
    for fcst in data:
        stat_dict = {
            'Fcst Hour': fcst,
            'AERONET_AOD_Mean': np.nanmean(data[fcst]['aod_550']),
            'AERONET_AOD_Std': np.nanstd(data[fcst]['aod_550']),
            'AERONET_Angstrom_Mean': np.nanmean(data[fcst]['angstrom']),
            'AERONET_Angstrom_Std': np.nanstd(data[fcst]['angstrom']),
        }
        
        for col in data[fcst].columns:
            if 'aod' in col and col not in ['aod_550', 'climatology_aod']:
                prefix = col.replace('_aod', '').upper()
                ang_col = col.replace('aod', 'ang') 
                
                stat_dict[f'{prefix}_AOD_Mean'] = np.nanmean(data[fcst][col])
                stat_dict[f'{prefix}_AOD_Std'] = np.nanstd(data[fcst][col])
                stat_dict[f'{prefix}_AOD_RMSE'] = np.sqrt(((data[fcst][col] - data[fcst]['aod_550'])**2).mean())
                stat_dict[f'{prefix}_AOD_Correlation'] = data[fcst][col].corr(data[fcst]['aod_550'])
                
                if ang_col in data[fcst].columns:
                    stat_dict[f'{prefix}_Angstrom_Mean'] = np.nanmean(data[fcst][ang_col])
                    stat_dict[f'{prefix}_Angstrom_Std'] = np.nanstd(data[fcst][ang_col])
                    stat_dict[f'{prefix}_Angstrom_RMSE'] = np.sqrt(((data[fcst][ang_col] - data[fcst]['angstrom'])**2).mean())
                    stat_dict[f'{prefix}_Angstrom_Correlation'] = data[fcst][ang_col].corr(data[fcst]['angstrom'])

                if 'climatology_aod' in data[fcst].columns:
                    col_anom = data[fcst][col] - data[fcst]['climatology_aod']
                    obs_anom = data[fcst]['aod_550'] - data[fcst]['climatology_aod']
                    stat_dict[f'{prefix}_AOD_ACC'] = col_anom.corr(obs_anom)
        if 'model_aod' in data[fcst].columns and 'analysis_aod' in data[fcst].columns:
            stat_dict['MODEL_ANALYSIS_AOD_RMSE'] = np.sqrt(((data[fcst]['model_aod'] - data[fcst]['analysis_aod'])**2).mean())
            stat_dict['MODEL_ANALYSIS_AOD_Correlation'] = data[fcst]['model_aod'].corr(data[fcst]['analysis_aod'])

        stats.append(stat_dict)
    return stats





#####################################################################################
# Main function
#####################################################################################
def main(start_date, end_date, base_path, experiment_name, output_dir="./aeronet_model_comparison/", min_points=30, ts_freq='none'):
    start_time = time.time()
    start_day = datetime.strftime(start_date, '%Y%m%d')
    start_hour = datetime.strftime(start_date, '%H')
    n_procs = max(1, mp.cpu_count() - 1)
    model_path_template = config['paths']['model1']
    aeronet_dir_base = config['paths']['observations']
    analysis_path_template = config['paths']['model2']
    file_comparison = 'aeronet' # file_comparison used for output file names

    
    print(f" Locating AERONET daily files between {start_date.strftime('%Y-%m-%d')} and {end_date.strftime('%Y-%m-%d')}...")
    
    files = []
    current_date = start_date.replace(hour=0, minute=0, second=0)

    num_analysis_files = file_types.count('analysis')
    analysis_counter = 1

    for i in range(len(file_paths)):
        if file_types[i] =='aeronet':
    
            ############################
            ### Processing aeronet files
            ############################
        
            while current_date <= end_date:
                YYYY = current_date.strftime('%Y') # setting date variables
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
        
            if all_times == True:
                full_time_range = pd.date_range(start=start_date, end=end_date, freq = 'h')
                station_info = master_df[['station','lats','lons']].drop_duplicates()
                idx = pd.MultiIndex.from_product([station_info['station'], full_time_range],
                                                 names = ['station', 'obs hour'])
                complete_grid = pd.DataFrame(index=idx).reset_index()
                complete_grid = pd.merge(complete_grid, station_info, on='station', how='left')
                master_df = pd.merge(complete_grid, master_df[['station', 'obs hour', 'aod_550', 'angstrom']], on=['station','obs hour'], how ='left')
                
            master_df['obs date'] = master_df['obs hour'].dt.floor('D')
            print(f"\nTotal valid AERONET hours across all stations: {len(master_df)}")
            print(f"\nSampling for {experiment_name}...")
            master_df['model_aod'] = np.nan
            master_df['model_ang'] = np.nan
            master_df['climatology_aod'] = np.nan
            master_df['initialization time'] = pd.NaT
            master_df['initialization hour'] = '00'
            master_df['fcst length'] = pd.Timedelta(hours=0)
        
            grouped = [(date, group) for date, group in master_df.groupby('obs date')]
        
        elif file_types[i] == 'forecast':
            #########################
            ### Processing model data
            #########################
        
            if "forecast" in config['evaluation']:
                file_comparison = file_comparison + '_fcst'
                process_func = partial(processmodel, model_path_template=model_path_template)
                processed_dfs = []
            
                with mp.Pool(n_procs) as pool:
                    print('This may take a while. Walk around and get a coffee while you wait.')
                    for i, result_df in enumerate(pool.imap_unordered(process_func, grouped)):
                        processed_dfs.append(result_df)
                        if (i+1) % 10 == 0 or (i+1) == len(grouped):
                            print(f"Processed {i+1}/{len(grouped)} model days...")
        
                valid_dfs = [df for df in processed_dfs if df is not None and not df.empty]
        
                if valid_dfs:
                    model_df = pd.concat(valid_dfs, ignore_index=True)
                else:
                    model_df = pd.DataFrame()
                final_df = model_df.copy()
            else:
                model_df = master_df.copy()
                final_df = master_df.copy()
                
        elif file_types[i] == 'analysis':
            ############################
            ### Processing analysis data
            ############################
        
        
            if "analysis" in config['evaluation']:
                
                file_comparison = file_comparison + '_analysis'
                current_analysis_path = file_paths[i]
                analysis_process_func = partial(processanalysis, analysis_path_template=current_analysis_path)
                analysis_processed_dfs = []
        
        
                with mp.Pool(n_procs) as pool:
                    for i, result_df in enumerate(pool.imap_unordered(analysis_process_func, grouped)):
                        analysis_processed_dfs.append(result_df)
                        if (i+1) % 10 == 0 or (i+1) == len(grouped):
                            print(f"Processed {i+1}/{len(grouped)} analysis days...")
        
                # cleaning data and combining data frames
                analysis_df = pd.concat(analysis_processed_dfs, ignore_index = True)
                if num_analysis_files > 1:
                    col_suffix = f"_{analysis_counter}"
                    analysis_counter += 1
                else:
                    col_suffix = ""
                analysis_subset = analysis_df[['station', 'obs hour','analysis_aod', 'analysis_ang']].rename(columns = {'analysis_aod':f'analysis_aod{col_suffix}','analysis_ang':f'analysis_ang{col_suffix}'})
                
                if 'final_df' not in locals():
                    final_df = master_df.copy()
                    
                    
                try:
                    final_df = pd.merge(model_df, analysis_subset, on = ['station', 'obs hour'], how = 'inner')
                except Exception as e:
                    final_df = pd.merge(master_df, analysis_subset, on = ['station', 'obs hour'], how = 'inner')
        
            valid_counts = final_df.dropna(subset=['aod_550'])['station'].value_counts()
            valid_stations = valid_counts[valid_counts >= min_points].index
            final_df = final_df[final_df['station'].isin(valid_stations)].copy()


    # creating netCDF file for analysis
    # file primarily used with plot_nc_summaryfigs.py to plot hourly forecast RMSE and BIAS
    
    nc_df = final_df.copy()
    nc_df['fcst_length'] = nc_df['fcst length'].dt.total_seconds()/ 3600
    nc_df = nc_df.rename(columns={'initialization time':'init_time', 'obs hour':'obs_hour'})
    base_cols = ['station', 'init_time', 'fcst_length', 'obs_hour', 'lats', 'lons', 'aod_550', 'angstrom']

    extra_cols = [col for col in nc_df.columns if 'model' in col or 'analysis' in col]
    nc_cols = base_cols + extra_cols
    
    # nc_cols = ['station', 'init_time', 'fcst_length', 'obs_hour', 'lats', 'lons', 'aod_550', 'angstrom', 'model_aod', 'model_ang', 'analysis_aod', 'analysis_ang']

    # print(nc_df.columns)
    nc_df = nc_df[nc_cols]
    nc_df = nc_df.groupby(['fcst_length', 'station','init_time']).first()
    ds = nc_df.to_xarray()
    date_str = f"{start_date.strftime('%Y%m%d')}_{end_date.strftime('%Y%m%d')}"
    nc_output_file = os.path.join(output_dir , f"{experiment_name}_{file_comparison}_measurements_{date_str}.nc")
    ds.to_netcdf(nc_output_file)
    print(f"NetCDF saved to: {nc_output_file}\n")


    #####################
    ### Calculating Stats
    #####################


    date_str = f"{start_date.strftime('%Y%m%d')}_{end_date.strftime('%Y%m%d')}"
    for hour in np.arange(0,24,int(config['time_step'][:-1])):
        houred_df = final_df.loc[final_df['initialization hour'] == f'{hour:02d}']
        
        for i in range(len(eval_bins)-1):
            fcst_hr = 'fcst_hr' + f'{eval_bins[i]:02d}'
            upper_bound = int(eval_bins[i] + int(config['eval_bins']))
            lower_bound = int(eval_bins[i])
            binned_eval[fcst_hr] = houred_df.loc[(houred_df['fcst length'] >= timedelta(hours=lower_bound)) & (houred_df['fcst length'] <  timedelta(hours=upper_bound))]

        # Get stats for each forecast hour and save *{x}hourly_{init_time}_stats*.csv
        hourly_stats = displayStats(binned_eval)
        if hourly_stats:
            hourly_stats_df = pd.DataFrame(hourly_stats)
            output_file = os.path.join(output_dir, f"{experiment_name}_{file_comparison}_{int(config['eval_bins']):02d}hourly_{hour:02d}z_stats_{date_str}.csv")
            hourly_stats_df.to_csv(output_file, index =False)
            print('\nHourly stats saved to', output_file, '\n')


    # get stats for each station and save *comparison_stats*.csv
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
            'ANALYSIS_Angstrom_Correlation': grp['analysis_ang'].corr(grp['angstrom']),


        })



    if stats:
        out_df = pd.DataFrame(stats) 
        output_file = os.path.join(output_dir , 
                                   f"{experiment_name}_{file_comparison}_comparison_stats_{ts_freq}_{date_str}.csv")
        out_df.to_csv(output_file, index=False)
        print(f"Stats saved to {output_file}")
    else:
        print(f"No stations met the minimum threshold of {min_points} valid matched points to calculate statistics.")

    # create a *timeseries*.csv for each station based on observation hour 
    if ts_freq != 'none' and not final_df.empty:
        print(f"\nGenerating {ts_freq} timeseries data...")
        
        cols_to_keep = ['obs hour', 'initialization time', 'station', 'lats', 'lons', 'aod_550', 'model_aod', 'analysis_aod', 
                        'angstrom', 'model_ang', 'analysis_ang']
        analysis_cols = [col for col in final_df.columns if 'analysis' in col]
        cols_to_keep.extend(analysis_cols)
        
        ts_df = final_df[cols_to_keep].copy()
        ts_df.rename(columns={'obs hour': 'time', 'aod_550': 'aeronet_aod', 'angstrom': 'aeronet_angstrom'}, inplace=True)
        ts_df.insert(0, 'experiment', experiment_name)
        
        if ts_freq == 'daily':
            ts_df['time'] = ts_df['time'].dt.floor('D')
            ts_df = ts_df.groupby(['experiment', 'station', 'time', 'lats', 'lons']).mean(numeric_only=True).reset_index()
        elif ts_freq == 'monthly':
            ts_df['time'] = ts_df['time'].dt.to_period('M').dt.to_timestamp()
            ts_df = ts_df.groupby(['experiment', 'station', 'time', 'lats', 'lons']).mean(numeric_only=True).reset_index()
            
        ts_output_file = os.path.join(output_dir , 
                                      f"{experiment_name}_{file_comparison}_timeseries_{ts_freq}_{date_str}.csv")
        
        ts_df.to_csv(ts_output_file, index=False)
        print(f"Time series data saved to {ts_output_file}")

    print(f"Total time: {(time.time() - start_time)/60:.1f} minutes!")

#####################################################################################
# Running main
#####################################################################################
if __name__ == "__main__":
    # to change settings for the date, temporal frequency, and what is being evaluatated, edit config.yaml
    # collecting settings from config.yaml
    with open('config.yaml', 'r') as file:
        config = yaml.safe_load(file)
        
    start_date = datetime.strptime(config['dates']['start'], '%Y-%m-%d %H:%M:%S') 
    end_date = datetime.strptime(config['dates']['end'], '%Y-%m-%d %H:%M:%S') 
    
    base_path = config['paths']['base']
    experiment_name = config['paths']['experiment']
    output_path = config['paths']['output']

    all_times = config['save_missing_obs']

    aeronet_path = config['paths']['observations']
    forecast_path = config['paths']['model1']
    FPanalysis_path = config['paths']['model2']

    file_paths = [aeronet_path, forecast_path, FPanalysis_path]
    file_types = config['evaluation']
    
    print(f'Settings from {file}') # When running multiple experiments, it may be useful to create multiple config files 
    print(f'Data collected from:\n{aeronet_path}\n{forecast_path}\n{FPanalysis_path}')

    output_dir = os.path.join(base_path, output_path, experiment_name)
    print(base_path)
    print('Output Dir:', output_dir)
    min_points = config['min_station_obs']   # of required hourly data points to be included in high-level stats
    ts_freq = config['temporal_res']
    evaluation = config['evaluation']
    time_step = config['time_step']
  
    col_names = ['obs hour','aod_550','angstrom','lats','lons','obs date','model_aod','model_ang',
                 'analysis_aod','analysis_ang','initialization time' ,'fcst length']
    
    time_period =  datetime.strptime(config['dates']['end'], '%Y-%m-%d %H:%M:%S') - datetime.strptime(config['dates']['start'],'%Y-%m-%d %H:%M:%S')
    time_period = int(time_period.total_seconds())/3600
    if timedelta(hours = time_period) < timedelta(hours=int(config['fcst_period'])):
        fcst_period = int(time_period)
    else:
        fcst_period = int(config['fcst_period'])


    eval_bins = np.arange(0,fcst_period +1, int(config['eval_bins']))

    binned_eval = {}

    main(
        start_date=start_date,
        end_date=end_date,
        base_path=base_path,
        experiment_name=experiment_name,
        output_dir=output_dir, 
        min_points=min_points, 
        ts_freq=ts_freq
    )
