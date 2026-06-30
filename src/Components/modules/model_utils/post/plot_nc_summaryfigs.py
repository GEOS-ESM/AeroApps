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

warnings.filterwarnings('ignore')

# --- GLOBALLY INCREASE BASE FONT SIZE ---
plt.rcParams.update({'font.size': 16})


#####################################################################################
# Plot Points Map
#####################################################################################
def plot_4panel_map(df, var_name, suptitle, bias_col1, bias_col2, corr_col1, corr_col2, output_path, var, vmax_val):
    """Generates a 4-panel map for Bias and Correlation comparing GEOS-CAM and GEOS-FP."""
    fig, axs = plt.subplots(2,2,figsize = (22,14), subplot_kw = {'projection': ccrs.PlateCarree()})
    axs = axs.flatten()
    #print(axs)
    #fig = plt.figure(figsize=(22, 14)) # Slightly larger figure to accommodate bigger fonts
    
    # Define projections and axes
    #axes = [plt.subplot(2, 2, i, projection=ccrs.PlateCarree(), height_ratios=1.5) for i in range(1, 5)]
    if var == 'BIAS':
        panel_labels = [
            '(a) GEOS-CAM Bias', '(b) GEOS-FP Bias', 
            '(c) GEOS-CAM Correlation', '(d) GEOS-FP Correlation'
        ]
        clabel = 'BIAS'
        vmins = [-vmax_val, -vmax_val, 0.0, 0.0]
        cmaps = [plt.cm.RdBu_r, plt.cm.RdBu_r, plt.cm.viridis, plt.cm.viridis]
    elif var == 'RMSE':
        panel_labels = [
            '(a) GEOS-CAM RMSE', '(b) GEOS-FP RMSE', 
            '(c) GEOS-CAM Correlation', '(d) GEOS-FP Correlation'
        ]
        clabel = 'RMSE'
        vmins = [0.0, 0.0, 0.0, 0.0]
        cmaps = [plt.cm.Reds, plt.cm.Reds, plt.cm.viridis, plt.cm.viridis]
    
    # Cap extremes for bias scaling (95th percentile)

    valid_data = df[[bias_col1, bias_col2]].dropna().values

    
    # if valid_data.size >0:
    #     max_bias_plot = np.percentile(np.abs(df[[bias_col1, bias_col2]].dropna().values), 95)
    # else:
    #     max_bias_plot = 0.0

    
    vmaxs = [vmax_val, vmax_val, 1.0, 1.0] # Correlation is capped at 1.0
    cols = [bias_col1, bias_col2, corr_col1, corr_col2]
    
    scatters = []
    
    for i in range(len(axs)):
        # print(df[cols[i]])
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
        
        # Increase Panel Title Font Size
        axs[i].set_title(panel_labels[i], fontsize=22, pad=12)
        axs[i].gridlines(draw_labels=False, alpha=0.2)
        
    plt.subplots_adjust(wspace=0.05, hspace=0.05)
    
    # Add Colorbars
    pos_top = axs[1].get_position()
    pos_bot = axs[3].get_position()
    
    # Bias Colorbar
    cbar_ax1 = fig.add_axes([pos_top.x1 + 0.02, pos_top.y0, 0.02, pos_top.height])
    cbar1 = plt.colorbar(scatters[1], cax=cbar_ax1)
    cbar1.set_label(f'{var_name} {clabel}', fontsize=20)
    cbar1.ax.tick_params(labelsize=16) # Increase tick number size
    
    # Correlation Colorbar
    cbar_ax2 = fig.add_axes([pos_bot.x1 + 0.02, pos_bot.y0, 0.02, pos_bot.height])
    cbar2 = plt.colorbar(scatters[3], cax=cbar_ax2)
    cbar2.set_label(f'{var_name} Correlation (r)', fontsize=20)
    cbar2.ax.tick_params(labelsize=16) # Increase tick number size
    
    # Increase Main Title Font Size
    fig.suptitle(f'{suptitle}\n{start_date} to {end_date}\nForecast Hour {int(fcst_hour):02d}', 
                 fontsize=26, fontweight='bold', y=0.96)
                 
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
    
    image_dir = "/discover/nobackup/caturne4/AeroApps/src/Components/modules/model_utils/post/comparison_figures/summary_map/"
    plot_type = ['RMSE', 'BIAS']
    
    for i in plot_type:
        if 'RMSE' == i:
            var = 'RMSE'
        elif 'BIAS' == i:
            var = 'BIAS'
        img_files = []
        for hr in range(0,73):
            file_name = f"aod_{var}_comparison_20260601_20260622_hr{hr:02d}_4panel_{init_hour:02d}z.png"
            # print(file_name)
            dates = file_name[20:37]
            # print(dates)
            img_files.append(image_dir +file_name)
            # search_pattern = os.path.join(image_folder, file_name)
            # print(search_pattern)
            # image_paths = sorted(glob.glob(img_files))
    
    
        images = [Image.open(img).convert("RGB") for img in img_files]
        
        # print(img_files)
        if images:
    
            images[0].save(
                f"./comparison_figures/summary_map/aod_{var}_{dates}_{init_hour}z.gif",     
                save_all=True,              
                append_images=images[1:],   
                duration=500,               
                loop=0                      
            )
            print(f"{var}_{dates}.gif created successfully with {len(images)} frames!")
        else:
            print(f"No images found matching the pattern:\n{search_pattern}")



#####################################################################################
# Main function
#####################################################################################
def main(ds, output_dir, fcst_hr, global_max_bias, global_max_rmse):
    os.makedirs(output_dir, exist_ok=True)


    
    ds = ds.sel(fcst_length = fcst_hr, method = 'nearest')
   
    
    
    # Filter for stations with at least 2000 valid data points
    original_len = len(ds.aod_550)

    df = pd.DataFrame({'Stations': ds['station'].values})

    # print(df)
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
    
    # df['MODEL_Angstrom_Bias'] = stn_model
    # df['ANALYSIS_Angstrom_Bias'] = ds['analysis_ang'] - ds['angstrom']

    df['MODEL_AOD_Correlation'] = xr.corr(ds['model_aod'], ds['aod_550'], dim = 'init_time')
    df['ANALYSIS_AOD_Correlation'] = xr.corr(ds['analysis_aod'], ds['aod_550'], dim = 'init_time')
    
    
    
    # 1. AOD 4-Panel Bias Plot
    plot_4panel_map(df, 
                    var_name="AOD",
                    suptitle="AOD at 550 nm Evaluated Against AERONET",
                    bias_col1='MODEL_AOD_Bias', bias_col2='ANALYSIS_AOD_Bias',
                    corr_col1='MODEL_AOD_Correlation', corr_col2='ANALYSIS_AOD_Correlation',
                    output_path=os.path.join(output_dir, f"aod_BIAS_comparison_{date_range}_hr{int(fcst_hour):02d}_4panel_{init_hour:02d}z.png"), 
                    var = 'BIAS', vmax_val = global_max_bias)
    #2. AOD 4-Panel RMSE Plot
    plot_4panel_map(df, 
                    var_name="AOD",
                    suptitle="AOD at 550 nm Evaluated Against AERONET",
                    bias_col1='MODEL_AOD_RMSE', bias_col2='ANALYSIS_AOD_RMSE',
                    corr_col1='MODEL_AOD_Correlation', corr_col2='ANALYSIS_AOD_Correlation',
                    output_path=os.path.join(output_dir, f"aod_RMSE_comparison_{date_range}_hr{int(fcst_hour):02d}_4panel_{init_hour:02d}z.png"),
                    var = 'RMSE', vmax_val = global_max_rmse)
                    
    # # 3. Angstrom 4-Panel Plot
    # plot_4panel_map(df, 
    #                 var_name="Angstrom Exponent",
    #                 suptitle="Angstrom Exponent Evaluated Against AERONET",
    #                 bias_col1='MODEL_Angstrom_Bias', bias_col2='ANALYSIS_Angstrom_Bias',
    #                 corr_col1='MODEL_Angstrom_Correlation', corr_col2='ANALYSIS_Angstrom_Correlation',
    #                 output_path=os.path.join(output_dir, "angstrom_comparison_4panel.png"), 
    #                 var = 'BIAS')
                    
    # # 4. Data Points Map
    # plot_points_map(df, output_path=os.path.join(output_dir, "data_points_coverage.png"))


#####################################################################################
# Running Main 
#####################################################################################
if __name__ == "__main__":

    file_dir = '/discover/nobackup/caturne4/geos_cam_eval/comparison/'
    n_points_min = 50
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', type=str, default=os.path.join(file_dir, '/discover/nobackup/caturne4/geos_cam_eval/comparison/geos_cam_eval_aeronet_fcst_analysis_measurements_20260601_20260622.nc'),
                        help='Path to the summary CSV file')
    parser.add_argument('--output-dir', type=str, default="./comparison_figures/summary_map/",
                        help='Directory to save output figures')
    args = parser.parse_args()
    date_range = args.input[-20:-3]
    start_date = date_range[:8]
    end_date = date_range[9:]
    # print(start_date, end_date)

    ds = xr.open_dataset(args.input)

    ds = ds.where((ds.lons >-132) & (ds.lons < -60))
    ds = ds.where((ds.lats >15) & (ds.lats < 51))

    init_hour = 12
    time_mask = ds.init_time.dt.hour == init_hour
    ds = ds.isel(init_time=time_mask)
    all_model_bias = (ds['model_aod'] - ds['aod_550']).mean(dim='init_time')
    all_analysis_bias = (ds['analysis_aod'] - ds['aod_550']).mean(dim='init_time')
    
    all_model_rmse = np.sqrt(((ds['model_aod'] - ds['aod_550'])**2).mean(dim='init_time'))
    all_analysis_rmse = np.sqrt(((ds['analysis_aod'] - ds['aod_550'])**2).mean(dim='init_time'))

    global_max_bias = np.nanpercentile(np.abs(xr.concat([all_model_bias, all_analysis_bias], dim='station')), 95)
    global_max_rmse = np.nanpercentile(np.abs(xr.concat([all_model_rmse, all_analysis_rmse], dim='station')), 95)
    


    
    fcst_hour = 72
    for i in range(0,fcst_hour+1):
        fcst_hour = i
        main(ds, args.output_dir, fcst_hour, global_max_bias, global_max_rmse)

    save_gif()
