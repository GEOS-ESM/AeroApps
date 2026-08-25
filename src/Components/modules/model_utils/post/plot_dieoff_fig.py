import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import argparse
import os
from datetime import datetime as dt
from datetime import timedelta
import xarray as xr


def dieoff(df, var):

    if df['ANALYSIS_AOD_RMSE'].isna().all():
        print('WARNING: ANALYSIS_AOD_RMSE is empty')
    if df['MODEL_AOD_RMSE'].isna().all():
        print('WARNING: MODEL_AOD_RMSE is empty')
    if df['ANALYSIS_AOD_Correlation'].isna().all():
        print('WARNING: ANALYSIS_AOD_Correlation is empty')
    if df['MODEL_AOD_Correlation'].isna().all():
        print('WARNING: MODEL_AOD_Correlation is empty')
    
    if var == 'RMSE':
        analysis_line = df['ANALYSIS_AOD_RMSE'] 
        model_line    = df['MODEL_AOD_RMSE']
        ytext = 'AOD RMSE'

    elif var == 'CORR':
        model_line    = df['MODEL_AOD_Correlation']
        analysis_line = df['ANALYSIS_AOD_Correlation']
        ytext = 'AOD Correlation'

    elif var == 'ACC':
        model_line    = df['MODEL_AOD_ACC']
        analysis_line = df['ANALYSIS_AOD_ACC']
        ytext = 'AOD Anomoly Correlation Coefficient'
        
    
    fcst_hr = df['Fcst Hour']
    hour1 = int(fcst_hr[0][7:11])
    hour2 = int(fcst_hr[1][7:11])
    time_step = hour2-hour1
    last_hr = int(fcst_hr[len(fcst_hr)-1][7:11])
    fcst_hours = np.arange(0,last_hr+1, time_step)
    # print(fcst_hours)

    fig = plt.figure(figsize = (7,4))
    plt.plot(fcst_hours, analysis_line, label = 'BASELINE')
    plt.plot(fcst_hours, model_line, label = 'HIGHRES')
    plt.legend()
    plt.xlabel('Forecast Hour', fontsize=16)
    plt.ylabel(ytext, fontsize=16)
    plt.grid(alpha=0.3)
 


def main():

    parser = argparse.ArgumentParser()
    base_dir = '/discover/nobackup/caturne4/geos_cam_eval/comparison/'
    #default_file = 'geos_cam_eval_+aeronet_fcst_6hourly_stats_20260601_20260615.csv'
    parser.add_argument('--corr_file', type=str, default ='geos_cam_eval_aeronet_fcst_analysis_6hourly_00z_stats_20260601_20260615.csv', help='Path to input file with correlations')

    args = parser.parse_args()
    file_path = os.path.join(base_dir, args.corr_file)
    time_index = file_path.find('hourly')
    time_scale = file_path[time_index-1:time_index+6]
    date_index = file_path.find('stats_')
    dates = file_path[date_index+6:date_index+23]
    init_hour_index = file_path.find('z')
    init_hour = file_path[init_hour_index-2:init_hour_index]
    
    # print(dates)

    if not os.path.exists(file_path):
        
        print(f"Warning: File '{file_path}' not found.")
        return None

    plots = ['RMSE', 'CORR']

    df = pd.read_csv(file_path)
    df['climate_aod'] = climate_mean
    # print(df['Fcst Hour'])
    # print(df.columns)
    for var in plots:

        dieoff(df, var)
        plt.title(f'{init_hour}Z AOD {var} Dieoff', fontsize = 20)
        img_name = f'./comparison_figures/test_aod_{var}_{time_scale}_{init_hour}z_{dates}_dieoff.png'
        plt.savefig(img_name)
        print(f'{img_name} saved successfully')

if __name__ == '__main__':
    main()

