#!/usr/bin/env python3
import numpy as np
import os
import sys
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature


def plot_4_panel_comparison(lon, lat, model_aod_mean_baseline, obs_aod_mean, model_aod_mean, EXPID,EXPID_baseline,sat_name):
    """
    Plots a 4-panel figure comparing:
    Panel A: Baseline Observed AOD (obs_aod_mean_baseline)
    Panel B: Current Observed AOD (obs_aod_mean)
    Panel C: Modeled AOD (model_aod_mean)
    Panel D: Closeness metric = |model AOD - obs AOD| - |baseline - obs AOD|

    Parameters:
        lon, lat (2D or 1D arrays): Longitude and Latitude data for the global grid
        obs_aod_mean_baseline (array): Baseline observational data (AOD)
        obs_aod_mean (array): Experimental observational data (AOD)
        model_aod_mean (array): AOD data from the simulation model
        EXPID (str): Experiment ID used for figure title
    """
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    
    # Compute Panel D: Closeness metric
    closeness = np.abs(model_aod_mean - obs_aod_mean) - np.abs(model_aod_mean_baseline - obs_aod_mean)
    
    # Create the figure with a 2x2 grid of panels
    fig, axes = plt.subplots(2, 2, figsize=(10, 6.5), subplot_kw={'projection': ccrs.PlateCarree()})
    
    # Common color settings for AOD plots
    vmin_aod = 0
    vmax_aod = 1
    vmin_closeness = -0.25
    vmax_closeness = 0.25
    cmap_aod = 'viridis'
    cmap_closeness = 'RdBu_r'
    
    # --- Panel A: Baseline Observed AOD ---
    im_a = axes[0, 0].pcolormesh(lon, lat, model_aod_mean_baseline, cmap=cmap_aod, vmin=vmin_aod, vmax=vmax_aod, transform=ccrs.PlateCarree())
    axes[0, 0].set_title(EXPID_baseline, fontsize=12)
    axes[0, 0].add_feature(cfeature.COASTLINE, edgecolor='black', linewidth=0.8)

    # --- Panel B: Current Observed AOD ---
    im_b = axes[0, 1].pcolormesh(lon, lat, obs_aod_mean, cmap=cmap_aod, vmin=vmin_aod, vmax=vmax_aod, transform=ccrs.PlateCarree())
    axes[0, 1].set_title(f"MODIS NNR ({sat_name})", fontsize=12)
    axes[0, 1].add_feature(cfeature.COASTLINE, edgecolor='black', linewidth=0.8)

    # --- Panel C: Modeled AOD ---
    im_c = axes[1, 0].pcolormesh(lon, lat, model_aod_mean, cmap=cmap_aod, vmin=vmin_aod, vmax=vmax_aod, transform=ccrs.PlateCarree())
    axes[1, 0].set_title(EXPID, fontsize=12)
    axes[1, 0].add_feature(cfeature.COASTLINE, edgecolor='black', linewidth=0.8)

    # --- Panel D: Closeness Analysis ---
    im_d = axes[1, 1].pcolormesh(lon, lat, closeness, cmap=cmap_closeness, vmin=vmin_closeness, vmax=vmax_closeness, transform=ccrs.PlateCarree())
    axes[1, 1].set_title("Closeness", fontsize=12)
    axes[1, 1].add_feature(cfeature.COASTLINE, edgecolor='black', linewidth=0.8)
    
    # Add colorbars for each row
    cax_c = fig.add_axes([0.1, 0.03, 0.3, 0.02])  # Colorbar for AOD rows
    cbar_c = fig.colorbar(im_c, cax=cax_c, orientation='horizontal', extend='both')
    cbar_c.set_label("AOD at 550 nm")

    cax_d = fig.add_axes([0.6, 0.03, 0.3, 0.02])  # Colorbar for closeness row
    cbar_d = fig.colorbar(im_d, cax=cax_d, orientation='horizontal', extend='both')
    cbar_d.set_label("|new-obs|-|baseline-obs|")

    # Main title for the figure
    #plt.suptitle(f"4-Panel Comparison: {EXPID}", fontsize=16, y=0.95)
    
    # Adjust layout
    plt.tight_layout(rect=[0, 0, 1, 0.93])
    plt.show()

def plot_GEOS_single_exp_compare(lon, lat, baseline_aod_mean,EXPID_baseline,target_aod_mean,EXPID,):
    fig, axes = plt.subplots(1, 3, figsize=(20, 6), subplot_kw={'projection': ccrs.PlateCarree()})
    
    # Common settings
    vmin_aod = 0  # Min value for AOD
    vmax_aod = 1  # Max value for AOD
    vmin_diff = -0.2  # Min value for difference
    vmax_diff = 0.2  # Max value for difference
    cmap_aod = 'viridis'  # Colormap for AOD
    cmap_diff = 'RdBu_r'  # Colormap for difference (bias)
    
    # Plot Observed AOD
    im1 = axes[0].pcolormesh(lon, lat, baseline_aod_mean, cmap=cmap_aod, vmin=vmin_aod, vmax=vmax_aod, transform=ccrs.PlateCarree())
    axes[0].set_title(EXPID_baseline)
    axes[0].add_feature(cfeature.COASTLINE, edgecolor='black', linewidth=1)
    #axes[0].add_feature(cfeature.BORDERS, edgecolor='black', linewidth=0.5)
    
    # Plot Modeled AOD
    im2 = axes[1].pcolormesh(lon, lat, target_aod_mean, cmap=cmap_aod, vmin=vmin_aod, vmax=vmax_aod, transform=ccrs.PlateCarree())
    axes[1].set_title(EXPID)
    axes[1].add_feature(cfeature.COASTLINE, edgecolor='black', linewidth=1)
    
    # Compute Bias 
    bias_aod = target_aod_mean - baseline_aod_mean
    # Plot Bias 
    im3 = axes[2].pcolormesh(lon, lat, bias_aod, cmap=cmap_diff, vmin=vmin_diff, vmax=vmax_diff, transform=ccrs.PlateCarree())
    axes[2].set_title("Difference (New - Baseline)")
    axes[2].add_feature(cfeature.COASTLINE, edgecolor='black', linewidth=1)
    
    
    # Add Colorbars for Each Plot
    cbar1 = fig.colorbar(im1, ax=axes[0], orientation='horizontal', pad=0.05, aspect=30,extend='both',shrink=0.7,)
    cbar1.set_label("AOD at 550 nm")
    
    cbar2 = fig.colorbar(im2, ax=axes[1], orientation='horizontal', pad=0.05, aspect=30,extend='both',shrink=0.7,)
    cbar2.set_label("AOD at 550 nm")
    
    cbar3 = fig.colorbar(im3, ax=axes[2], orientation='horizontal', pad=0.05, aspect=30,extend='both',shrink=0.7,)
    cbar3.set_label("Difference")
    plt.suptitle(f"TOTEXTTAU550", fontsize=14, y=0.85)
    
    # Display the plot
    plt.tight_layout()
    plt.show()


def read_model(EXPID,EXPDIR,yy,mm,sat,varn="TOTEXTTAU550"):
    varout = []
    lat, lon = None, None
    filen = []
    
    for m in mm:
        filen_ = f"{EXPDIR}/holding/inst2d_hwl_x/{yy:04d}{m:02d}/{EXPID}.inst2d_hwl_x.monthly.{sat}.{yy:04d}{m:02d}.nc4"
        if not os.path.exists(filen_):
            print(f"File not found: {filen_}. Skipping...")
            continue
        filen.append(filen_)
        print(filen)
    
    ds = xr.open_mfdataset(filen)
    varout = ds[varn].mean(dim='time')
    lat = ds['lat'].values
    lon = ds['lon'].values
        
    return varout, lon, lat

def read_sat(sat,yy,mm,varn="tau_"):
    varout = []
    lat, lon = None, None
    filen = []
    
    for m in mm:
        filen_ = f"/home/pcolarco/geos_aerosols/pcolarco/{sat}/Level3/Y{yy:04d}/M{m:02d}/nnr_003.{sat}_L3a.blend.monthly.{yy:04d}{m:02d}.nc4"
        if not os.path.exists(filen_):
            print(f"File not found: {filen_}. Skipping...")
            continue
        filen.append(filen_)
        print(filen)
    
    ds = xr.open_mfdataset(filen)
    varout = ds[varn].mean(dim='time')
    lat = ds['lat'].values
    lon = ds['lon'].values
        
    return varout, lon, lat
    

if __name__ == "__main__":
    baseline = "c180R_v11.8.0_develop"
    perturb  = "c180R_v11.8.0_newdust"
    sat      = "MOD04"
    yy       = 2019
    mm       = [1,2]
    basedir  = f"/home/pcolarco/geos_aerosols/pcolarco/{baseline}/"
    pertdir  = f"/home/pcolarco/geos_aerosols/pcolarco/{perturb}/"
    baseaod, lon, lat  = read_model(baseline,basedir,yy,mm,sat)
    pertaod, lon, lat  = read_model(perturb ,pertdir,yy,mm,sat)
    obsaod, lono, lato = read_sat(sat,yy,mm)
    #plot_GEOS_single_exp_compare(lon, lat, baseaod.values,baseline,pertaod.values,perturb,)
    obsaodi = obsaod.interp_like(baseaod,method="nearest",assume_sorted=True)
    #plot_GEOS_single_exp_compare(lon, lat, baseaod.values,baseline,np.squeeze(obsaodi.values),sat,)
    plot_4_panel_comparison(lon, lat, baseaod.values, np.squeeze(obsaodi.values), pertaod.values, perturb,baseline,sat)
