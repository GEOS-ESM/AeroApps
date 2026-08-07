'''
This code will plot a 6 panel figure (2 rows by 3 columns) between two models.
It will have bias on the left, rmse in the middle, and correlation on the right.
Each plot contains the average bias, rmse, and correlation for a 12 hour time chunk.
This is similar to plot_6panel_summaryfigs_vert.py
'''


import os
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import argparse
import warnings
import glob
from PIL import Image
import math

warnings.filterwarnings('ignore')

# --- GLOBALLY INCREASE BASE FONT SIZE ---
plt.rcParams.update({'font.size': 16})


#####################################################################################
# Plot Points Map
#######################z#############################################################
def plot_6panel_map(df, var_name, suptitle, bias_col1, bias_col2, rmse_col1, rmse_col2, corr_col1, corr_col2, output_path, var, bias_vmax_val, rmse_vmax_val):
    """Generates a 4-panel map for Bias and Correlation comparing GEOS-CAM and GEOS-FP."""
    fig, axs = plt.subplots(2,3,figsize = (27,10), subplot_kw = {'projection': ccrs.PlateCarree()})
    axs = axs.flatten()
    colorbar_buffer = 0.03
    bias_vmax_val = math.ceil((bias_vmax_val+colorbar_buffer)*100)/100
    bias_vmax_val = round(round(bias_vmax_val*20)/20,2)
    rmse_vmax_val = math.ceil((rmse_vmax_val+colorbar_buffer)*100)/100
    rmse_vmax_val = round(round(rmse_vmax_val*20)/20,2)
    
    
    panel_labels = [
        '(a) GEOS-CAM Bias', '(b) GEOS-CAM RMSE','(c) GEOS-CAM Correlation',
        '(d) GEOS-FP Bias', '(e) GEOS-FP RMSE', '(f) GEOS-FP Correlation'
    ]
    
    vmins = [-bias_vmax_val, 0.0,0.0, -bias_vmax_val, 0.0, 0.0]
    cmaps = [plt.cm.RdBu_r,plt.cm.Reds, plt.cm.viridis, plt.cm.RdBu_r, plt.cm.Reds, plt.cm.viridis]
    bias_clabel = 'BIAS'
    
    rmse_clabel = 'RMSE'


    valid_data = df[[bias_col1, bias_col2, rmse_col1, rmse_col2]].dropna().values

    
    vmaxs = [bias_vmax_val, rmse_vmax_val, 1.0, bias_vmax_val, rmse_vmax_val, 1.0] # Correlation is capped at 1.0
    cols = [bias_col1, rmse_col1, corr_col1, bias_col2, rmse_col2, corr_col2]
    
    scatters = []
    
    for i in range(len(axs)):

        axs[i].add_feature(cfeature.COASTLINE, linewidth=0.5)
        axs[i].add_feature(cfeature.BORDERS, linestyle=':', linewidth=0.5)
        axs[i].add_feature(cfeature.LAND, facecolor='lightgray', alpha=0.3)
        axs[i].add_feature(cfeature.STATES, alpha = 0.3)
        axs[i].set_extent([-126,-65,24,50], crs = ccrs.PlateCarree())
        clean_df = df.dropna(subset = ['Lons', 'Lats', cols[i]])
        sc = axs[i].scatter(clean_df['Lons'], clean_df['Lats'], c=clean_df[cols[i]], 
                        cmap=cmaps[i], vmin=vmins[i], vmax=vmaxs[i],
                        s=80, edgecolor='k', linewidth=0.3, transform=ccrs.PlateCarree())
        scatters.append(sc)

        axs[i].set_title(panel_labels[i], fontsize=22, pad=12)
        axs[i].gridlines(draw_labels=False, alpha=0.2)
        
    plt.subplots_adjust(wspace=0.1, hspace=0.05)
    
    # Add Colorbars
    pos_top = axs[3].get_position()
    pos_mid = axs[4].get_position()
    pos_bot = axs[5].get_position()
    
    # Bias or RMSE Colorbar    left, bottom, width, height
    cbar_ax1 = fig.add_axes([pos_top.x0, pos_top.y0- 0.04, pos_top.width,0.03])
    cbar1 = plt.colorbar(scatters[3], cax=cbar_ax1, orientation = 'horizontal')

    cbar1.set_label(f'{var_name} {bias_clabel}', fontsize=20)
    cbar1.ax.tick_params(labelsize=16) # Increase tick number size
    
    # Correlation Colorbar
    cbar_ax2 = fig.add_axes([pos_mid.x0, pos_mid.y0 - 0.04, pos_mid.width, 0.03])
    cbar2 = plt.colorbar(scatters[4], cax=cbar_ax2, orientation = 'horizontal')
    cbar2.set_label(f'{var_name} {rmse_clabel}', fontsize=20)
    cbar2.ax.tick_params(labelsize=16) # Increase tick number size
    
    cbar_ax2 = fig.add_axes([pos_bot.x0, pos_bot.y0 - 0.04, pos_bot.width, 0.03])
    cbar2 = plt.colorbar(scatters[5], cax=cbar_ax2, orientation = 'horizontal')
    cbar2.set_label(f'{var_name} Correlation (r)', fontsize=20)
    cbar2.ax.tick_params(labelsize=16) # Increase tick number size
    
    if fcst_hour == 0:

        fig.suptitle(f'Forecast Hour {int(fcst_hour):02d}',fontsize=25, fontweight='bold', y=0.93)
    else:
        fig.suptitle(f'Forecast Hours {int(fcst_hour-12):02d} - {int(fcst_hour):02d}',fontsize=25, fontweight='bold', y=0.93)
        
                 
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: {output_path}")



#####################################################################################
# Plot Points Map
#####################################################################################
def plot_points_map(df, output_path):
    """Generates a single-panel map showing the number of matched data points."""
    fig = plt.figure(figsize=(14, 9))
    ax = plt.subplot(1, 1, 1, projection=ccrs.PlateCarree())
    
    ax.add_feature(cfeature.COASTLINE)
    ax.add_feature(cfeature.BORDERS, linestyle=':')
    ax.add_feature(cfeature.STATES, alpha = 0.3)
    ax.add_feature(cfeature.LAND, facecolor='lightgray', alpha=0.3)
    ax.set_extent([-132,-60,23,60])
    
    # Use log scale for point sizes to handle large variations
    sizes = 40 + 150 * (np.log10(df['N_Points']) - np.log10(df['N_Points'].min())) / \
            (np.log10(df['N_Points'].max()) - np.log10(df['N_Points'].min()))
            
    sc = ax.scatter(df['Lons'], df['Lats'], c=df['N_Points'], 
                    cmap='magma_r', s=sizes, edgecolor='k', linewidth=0.5, 
                    transform=ccrs.PlateCarree(), alpha=0.8)
                    
    cbar = plt.colorbar(sc, shrink=0.6, pad=0.04)
    cbar.set_label('Number of Matched Hourly Observations', fontsize=18)
    cbar.ax.tick_params(labelsize=16)
    
    # Increase Title Font Size
    ax.set_title(f'Valid AERONET vs HIGHRES MODEL Data Points\n({len(df)} stations)', 
                 fontsize=22, pad=15)
                 
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: {output_path}")

# 
def save_gif():
    
    image_dir = "/discover/nobackup/caturne4/AeroApps/src/Components/modules/model_utils/post/final_analysis/time_chunks/"
    
    
    img_var = 'BIAS'
    save_var = 'ALL'
    img_files = []
    for hr in range(0,fcst_hour+1,12):
        file_name = f"aod_{img_var}_comparison_{start_date}_{end_date}_hr{hr:02d}_6panel_{init_hour:02d}z.png"
        dates = file_name[20:37]
        img_files.append(image_dir +file_name)

    images = [Image.open(img).convert("RGB") for img in img_files]
    
    if images:

        images[0].save(os.path.join(image_dir,
            f"aod_6panel_{save_var}_{dates}_{init_hour:02d}z.gif"),     
            save_all=True,              
            append_images=images[1:],   
            duration=2500,               
            loop=0                      
        )
        print(f"{file_name} created successfully with {len(images)} frames!")
    else:
        print(f"No images found matching the pattern:\n{search_pattern}")


#####################################################################################
# Main function
#####################################################################################
def main(full_ds, output_dir, start_hour, fcst_hr, global_max_bias, global_max_rmse):
    os.makedirs(output_dir, exist_ok=True)
    if fcst_hr == 0:
        ds = full_ds.sel(fcst_length = fcst_hr, method = 'nearest')
    else:

        sliced_ds = full_ds.sel(fcst_length = slice(fcst_hr -12,fcst_hr))
        ds = sliced_ds.mean(dim = 'fcst_length', skipna = True)
   
    original_len = len(ds.aod_550)

    df = pd.DataFrame({'Stations': ds['station'].values})

    model_bias = ds['model_aod'] - ds['aod_550'] 
    stn_model_bias = model_bias.mean(dim = 'init_time')
    stn_model_rmse = np.sqrt((model_bias**2).mean(dim='init_time'))
    
    analysis_bias = ds['analysis_aod'] - ds['aod_550'] 
    stn_analysis_bias = analysis_bias.mean(dim = 'init_time')
    stn_analysis_rmse = np.sqrt((analysis_bias**2).mean(dim='init_time'))

    
    # Calculate Biases
    df['MODEL_AOD_Bias'] = stn_model_bias
    df['ANALYSIS_AOD_Bias'] = stn_analysis_bias
    
    df['MODEL_AOD_RMSE'] = stn_model_rmse
    df['ANALYSIS_AOD_RMSE'] = stn_analysis_rmse
    df['Lats'] = ds['lats'].mean(dim='init_time')
    df['Lons'] = ds['lons'].mean(dim='init_time')

    df['MODEL_AOD_Correlation'] = xr.corr(ds['model_aod'], ds['aod_550'], dim = 'init_time')
    df['ANALYSIS_AOD_Correlation'] = xr.corr(ds['analysis_aod'], ds['aod_550'], dim = 'init_time')
    
    plot_6panel_map(df, 
                    var_name="AOD",
                    suptitle="AOD at 550 nm Evaluated Against AERONET",
                    bias_col1='MODEL_AOD_Bias', bias_col2='ANALYSIS_AOD_Bias',
                    rmse_col1='MODEL_AOD_RMSE', rmse_col2='ANALYSIS_AOD_RMSE',
                    corr_col1='MODEL_AOD_Correlation', corr_col2='ANALYSIS_AOD_Correlation',
                    output_path=os.path.join(output_dir, f"aod_horiz_BIAS_comparison_{date_range}_hr{int(fcst_hour):02d}_6panel_{init_hour:02d}z.png"), 
                    var = 'BIAS', bias_vmax_val = global_max_bias, rmse_vmax_val = global_max_rmse)
    


#####################################################################################
# Running Main 
#####################################################################################
if __name__ == "__main__":

    file_dir = '/discover/nobackup/caturne4/geos_cam_eval/comparison/'
    n_points_min = 30
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', type=str, default=os.path.join(file_dir, '/discover/nobackup/caturne4/AeroApps/src/Components/modules/model_utils/post/final_analysis/mix/mix_aeronet_fcst_analysis_measurements_20260529_20260727.nc'),


                        help='Path to the summary CSV file')
    parser.add_argument('--output-dir', type=str, default="./time_chunks/",
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


    ds = xr.open_dataset(args.input)

    ds = ds.where((ds.lons >-132) & (ds.lons < -60))
    ds = ds.where((ds.lats >15) & (ds.lats < 51))

    for init_hour in [0, 12]:
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

        start_hour = 0
    
        for i in range(0,fcst_hour+1,12):
            
            fcst_hour = i
            main(filtered_ds, args.output_dir, start_hour, fcst_hour, global_max_bias, global_max_rmse)
            start_hour +=12
    
        save_gif()
