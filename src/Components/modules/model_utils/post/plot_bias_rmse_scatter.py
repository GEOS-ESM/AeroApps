import os
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
import argparse
import warnings
import glob
from PIL import Image
import math
import warnings
warnings.filterwarnings('ignore')
import matplotlib.cm as cm
import matplotlib.colors as mcolors
from datetime import datetime

def main(ds, output_dir, fcst_hr, var_max, var_min, color):
    os.makedirs(output_dir, exist_ok=True)    
    ds = ds.sel(fcst_length = fcst_hr, method = 'nearest')
    df = pd.DataFrame({'Stations': ds['station'].values})
 
    model_bias = ds['model_aod'] - ds['aod_550'] 
    stn_model_bias = model_bias.mean(dim = 'init_time')
    stn_model_rmse = np.sqrt((model_bias**2).mean(dim='init_time'))
    
    analysis_bias = ds['analysis_aod'] - ds['aod_550'] 
    stn_analysis_bias = analysis_bias.mean(dim = 'init_time')
    stn_analysis_rmse = np.sqrt((analysis_bias**2).mean(dim='init_time'))

    
    # Calculate Biases and RMSE
    df['MODEL_AOD_Bias'] = stn_model_bias
    df['ANALYSIS_AOD_Bias'] = stn_analysis_bias
    
    df['MODEL_AOD_RMSE'] = stn_model_rmse
    df['ANALYSIS_AOD_RMSE'] = stn_analysis_rmse
    df['Lats'] = ds['lats'].mean(dim='init_time')
    df['Lons'] = ds['lons'].mean(dim='init_time')

    df['MODEL_AOD_Correlation'] = xr.corr(ds['model_aod'], ds['aod_550'], dim = 'init_time')
    df['ANALYSIS_AOD_Correlation'] = xr.corr(ds['analysis_aod'], ds['aod_550'], dim = 'init_time')

    
    if var == 'BIAS':
        if var_max < np.nanmax(df['MODEL_AOD_Bias']):
            var_max = np.nanmax(df['MODEL_AOD_Bias'])
        if var_min > np.nanmin(df['MODEL_AOD_Bias']):
            var_min = np.nanmin(df['MODEL_AOD_Bias'])
            
        if var_max < np.nanmax(df['ANALYSIS_AOD_Bias']):
            var_max = np.nanmax(df['ANALYSIS_AOD_Bias'])
            # print(var_max)
        if var_min > np.nanmin(df['ANALYSIS_AOD_Bias']):
            var_min = np.nanmin(df['ANALYSIS_AOD_Bias'])
        plt.scatter(df['MODEL_AOD_Bias'], df['ANALYSIS_AOD_Bias'], color = color)
        x = df['MODEL_AOD_Bias']
        y = df['ANALYSIS_AOD_Bias']
    elif var == 'RMSE':
        if var_max < np.nanmax(df['MODEL_AOD_RMSE']):
            var_max = np.nanmax(df['MODEL_AOD_RMSE'])
        if var_min > np.nanmin(df['MODEL_AOD_RMSE']):
            var_min = np.nanmin(df['MODEL_AOD_RMSE'])
            
        if var_max < np.nanmax(df['ANALYSIS_AOD_RMSE']):
            var_max = np.nanmax(df['ANALYSIS_AOD_RMSE'])
        if var_min > np.nanmin(df['ANALYSIS_AOD_RMSE']):
            var_min = np.nanmin(df['ANALYSIS_AOD_RMSE'])

        ax.annotate('CAM Performs\nBetter', 
             xy=(0.1,0.75), 
             xytext=(0.1, 0.75),
             fontsize=15,
             color='black',
             bbox=dict(boxstyle='round,pad=0.5', fc='lightgrey', alpha=0.5))
        ax.annotate('FP Performs\nBetter', 
             xy=(0.75,0.05), 
             xytext=(0.75,0.05),
             fontsize=15,
             color='black',
             bbox=dict(boxstyle='round,pad=0.5', fc='lightgrey', alpha=0.1))
        
        ax.scatter(df['MODEL_AOD_RMSE'], df['ANALYSIS_AOD_RMSE'], color=color)
 
        x = df['MODEL_AOD_RMSE']
        y = df['ANALYSIS_AOD_RMSE']

    elif var == 'CORR':
        if var_max < np.nanmax(df['MODEL_AOD_Correlation']):
            var_max = np.nanmax(df['MODEL_AOD_Correlation'])
        if var_min > np.nanmin(df['MODEL_AOD_Correlation']):
            var_min = np.nanmin(df['MODEL_AOD_Correlation'])
            
        if var_max < np.nanmax(df['ANALYSIS_AOD_Correlation']):
            var_max = np.nanmax(df['ANALYSIS_AOD_Correlation'])
        if var_min > np.nanmin(df['ANALYSIS_AOD_Correlation']):
            var_min = np.nanmin(df['ANALYSIS_AOD_Correlation'])
        plt.scatter(df['MODEL_AOD_Correlation'], df['ANALYSIS_AOD_Correlation'], color = color)
        x = df['MODEL_AOD_Correlation']
        y = df['ANALYSIS_AOD_Correlation']
    
    return var_max, var_min, x, y

#####################################################################################
# Running Main 
#####################################################################################
if __name__ == "__main__":

    file_dir = '/discover/nobackup/caturne4/AeroApps/src/Components/modules/model_utils/post/final_analysis/'
    n_points_min = 30
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', type=str, default=os.path.join(file_dir, 'final_analysis_aeronet_forecast_analysis_forecast_analysis__measurements_20260715_20260727.nc'),
                        help='Path to the summary CSV file')
    parser.add_argument('--output-dir', type=str, default="./bias_rmse/",
                        help='Directory to save output figures')
    args = parser.parse_args()
    date_range = args.input[-20:-3]
    
    start_date = date_range[:8]
    start_year = start_date[:4]
    start_month = start_date[4:6]
    start_day = start_date[6:]
    
    end_date = date_range[9:]
    end_year = end_date[:4]
    end_month = end_date[4:6]
    end_day = end_date[6:]

    stations = ['ALL', 'GSFC', 'NEON_UNDE','Granite_Island','Chiwaukee_Prairie','NPU_Chicago_IL',
            'London-CDN','Egbert','Toronto', 'Dayton', 'East_Brunswick', 'Flushing','UMBC','Billerica']

    full_ds = xr.open_dataset(args.input)

    full_ds = full_ds.where((full_ds.lons >-132) & (full_ds.lons < -60))
    full_ds = full_ds.where((full_ds.lats >15) & (full_ds.lats < 51))

    portland_colors = ['#0c3383', '#0a88ba', '#f2d338', '#f28f38', '#d91e1e']
    portland_cmap = mcolors.LinearSegmentedColormap.from_list('portland', portland_colors)
    cmap = portland_cmap
    norm = mcolors.Normalize(vmin = 0, vmax = 72)

    
    for stn in stations:
        if stn != 'ALL':
            ds = full_ds.where(full_ds['station'] == stn)
            lat = ds['lats'].min()
            lon = ds['lons'].min()
        else:
            ds = full_ds
            
        for var in ['BIAS', 'RMSE']:
            var_max = -999
            var_min = 999
            if var == 'RMSE':
                fig, ax = plt.subplots(figsize = (12,9))
                ax.set_xlim(-0.05,var_max+0.05)
                ax.set_ylim(-0.05,var_max+0.05)

                ax.set_xlabel(f'GEOS-CAM', fontsize = 20)
                ax.set_ylabel(f'GEOS-FP', fontsize = 20)
                try:
                    ax.set_title(f'{stn}\nLat: {float(lat):.2f} Lon: {float(lon):.2f}', loc = 'right', fontsize = 16)
                    ax.set_title(f'GEOS-CAM vs GEOS-FP {var}\n{datetime.strptime(start_date, "%Y%m%d").date()} to {datetime.strptime(end_date, "%Y%m%d").date()}', loc = 'left', fontsize = 24)
                except:
                    ax.set_title(f'GEOS-CAM vs GEOS-FP {var} - {stn}\n{datetime.strptime(start_date, "%Y%m%d").date()} to {datetime.strptime(end_date, "%Y%m%d").date()}', fontsize = 24)
            elif var == 'BIAS':
                fig = plt.figure(figsize = (12,9))
                plt.xlabel(f'GEOS-CAM', fontsize = 20)
                plt.ylabel(f'GEOS-FP', fontsize = 20)
                plt.xlim(var_min-0.05,var_max+0.05)
                plt.ylim(var_min-0.05,var_max+0.05)
                try:
                    plt.title(f'{stn}\nLat: {float(lat):.2f} Lon: {float(lon):.2f}', loc = 'right', fontsize = 16)
                    plt.title(f'GEOS-CAM vs GEOS-FP {var}\n{datetime.strptime(start_date, "%Y%m%d").date()} to {datetime.strptime(end_date, "%Y%m%d").date()}', loc = 'left', fontsize = 24)
                except:
                    plt.title(f'GEOS-CAM vs GEOS-FP {var} - {stn}\n{datetime.strptime(start_date, "%Y%m%d").date()} to {datetime.strptime(end_date, "%Y%m%d").date()}', fontsize = 24)

            elif var == 'CORR':
                fig = plt.figure(figsize = (12,9))
                plt.xlabel(f'GEOS-CAM', fontsize = 20)
                plt.ylabel(f'GEOS-FP', fontsize = 20)
                plt.xlim(var_min-0.05,var_max+0.05)
                plt.ylim(var_min-0.05,var_max+0.05)
                try:
                    plt.title(f'{stn}\nLat: {float(lat):.2f} Lon: {float(lon):.2f}', loc = 'right', fontsize = 16)
                    plt.title(f'GEOS-CAM vs GEOS-FP {var}\n{datetime.strptime(start_date, "%Y%m%d").date()} to {datetime.strptime(end_date, "%Y%m%d").date()}', loc = 'left', fontsize = 24)
                except:
                    plt.title(f'GEOS-CAM vs GEOS-FP {var} - {stn}\n{datetime.strptime(start_date, "%Y%m%d").date()} to {datetime.strptime(end_date, "%Y%m%d").date()}', fontsize = 24)
    
            x_data = []
            y_data = []
            try: 
                plt.xticks(fontsize = 16)
                plt.yticks(fontsize = 16)
            except:
                ax.tick_params(axis='both', labelsize=16)
            for init_hour in [0,12]:

                time_mask = ds.init_time.dt.hour == init_hour
                filtered_ds = ds.isel(init_time=time_mask)
        
                all_model_bias = (filtered_ds['model_aod'] - filtered_ds['aod_550']).mean(dim='init_time')
                all_analysis_bias = (filtered_ds['analysis_aod'] - filtered_ds['aod_550']).mean(dim='init_time')
                
                all_model_rmse = np.sqrt(((filtered_ds['model_aod'] - filtered_ds['aod_550'])**2).mean(dim='init_time'))
                all_analysis_rmse = np.sqrt(((filtered_ds['analysis_aod'] - filtered_ds['aod_550'])**2).mean(dim='init_time'))
            
                global_max_bias = np.nanpercentile(np.abs(xr.concat([all_model_bias, all_analysis_bias], dim='station')), 95)
                global_max_rmse = np.nanpercentile(np.abs(xr.concat([all_model_rmse, all_analysis_rmse], dim='station')), 95)
            
        
                if init_hour == 0:
                    fcst_hour = 72
                elif init_hour == 12:
                    fcst_hour = 48
            
                for i in range(0,fcst_hour+1):
                    fcst_hour = i
                    color =cmap(norm(fcst_hour))
                    if stn!='ALL':
                        var_max, var_min, x_vals, y_vals = main(filtered_ds, args.output_dir, fcst_hour, var_max, var_min, color)
                    else:
                        var_max, var_min, x_vals, y_vals = main(filtered_ds, args.output_dir, fcst_hour, var_max, var_min, color)
                        if var == 'RMSE':
                            var_max = 1.25
                            var_min = -0.05
                        elif var =='BIAS':
                            var_max =0.55
                            var_min = -0.55
                    x_data.extend(np.ravel(x_vals))
                    y_data.extend(np.ravel(y_vals))
    
            if var == 'RMSE':
                ax.set_xlim(-0.05,var_max+0.05)
                ax.set_ylim(-0.05,var_max+0.05)
                valid_mask = ~np.isnan(x_data) & ~np.isnan(y_data)
                clean_x = np.array(x_data)[valid_mask]
                clean_y = np.array(y_data)[valid_mask]
        
                if len(clean_x) > 0:
                    slope, intercept = np.polyfit(clean_x, clean_y, 1)
                    line_x = np.array([var_min, var_max])
                    line_y = slope * line_x + intercept
                    plt.plot(line_x, line_y, color='orangered', linewidth=2, label=f'Regression (m={slope:.2f})')
                plt.plot([var_min-0.05, var_max+0.05], [var_min-0.05, var_max+0.05], 'k--', label='1:1 line')
                ax.legend(fontsize = 'large')
            elif var == 'BIAS':
                plt.xlim(var_min-0.05,var_max+0.05)
                plt.ylim(var_min-0.05,var_max+0.05)
                plt.plot([var_min-0.05, var_max+0.05], [var_min-0.05, var_max+0.05], 'k--', label='1:1 line')
                valid_mask = ~np.isnan(x_data) & ~np.isnan(y_data)
                clean_x = np.array(x_data)[valid_mask]
                clean_y = np.array(y_data)[valid_mask]
                if len(clean_x) > 0:
                    slope, intercept = np.polyfit(clean_x, clean_y, 1)
                    line_x = np.array([var_min, var_max])
                    line_y = slope * line_x + intercept
                plt.plot(line_x, line_y, color='orangered', linewidth=2, label=f'Regression (m={slope:.2f})')
                plt.legend(fontsize = 16)


            elif var == 'CORR':
                plt.xlim(var_min-0.05,var_max+0.05)
                plt.ylim(var_min-0.05,var_max+0.05)
                plt.plot([var_min-0.05, var_max+0.05], [var_min-0.05, var_max+0.05], 'k--', label='1:1 line')
                valid_mask = ~np.isnan(x_data) & ~np.isnan(y_data)
                clean_x = np.array(x_data)[valid_mask]
                clean_y = np.array(y_data)[valid_mask]
                if len(clean_x) > 0:
                    slope, intercept = np.polyfit(clean_x, clean_y, 1)
                    line_x = np.array([var_min, var_max])
                    line_y = slope * line_x + intercept
                plt.plot(line_x, line_y, color='orangered', linewidth=2, label=f'Regression (m={slope:.2f})')
                plt.legend(fontsize = 16)
    
            plt.grid(zorder = 0, alpha = 0.45)
            sm = cm.ScalarMappable(cmap=cmap, norm=norm)
            sm.set_array([])  
            cbar = fig.colorbar(sm, ax=plt.gca())
            cbar.set_label('Forecast Hour', fontsize=18)
            cbar.ax.tick_params(labelsize = 16)
            fig_path = os.path.join(args.output_dir, f'case_study_{stn}_{var}_{date_range}_scatter')
            plt.savefig(fig_path)
            print(f'made {stn} {var} plot')
        
