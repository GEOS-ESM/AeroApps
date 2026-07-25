#!/usr/bin/env python3
import argparse
import xarray as xr
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from cartopy import config
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from datetime import timezone
from astral import Observer
from astral.sun import sun
import os
os.environ['CARTOPY_USER_BACKGROUNDS'] = "/home/pcolarco/silo/python/"

import warnings
warnings.filterwarnings("ignore", category=UserWarning)

def main():
    parser = argparse.ArgumentParser(description='Plot pyroCb-cloud tracks from NetCDF file.')
    parser.add_argument('--event', type=str, default='Pas', 
                        help='Name of the event (default: Pas). Used to find Track_<event>.nc')
    args = parser.parse_args()
    
    event_name = args.event
   

    # Define custom file paths
    input_file = f'./Save/Track_{event_name}.nc'
    output_plot = f'./Plot/{event_name}_Plot_Refined.png'
    output_text = f'./Plot/{event_name}_Summary.txt'


    # 1. Load data and convert to DataFrame
    try:
        ds = xr.open_dataset(input_file)
        df = ds.to_dataframe()
    except FileNotFoundError:
        print(f"Error: Could not find file {input_file}")
        return

    # 2. FILTER OUT THE -1 BACKGROUND CELL
    df_clean = df[df['cell'] != -1].copy()

    # 3. IDENTIFY VALID TRACKS (> 2 HOURS & > 2 DEGREES)
    valid_cells = []
    for cell_id, group in df_clean.groupby('cell'):
        # Time filter
        duration = group['time'].max() - group['time'].min()
        if isinstance(duration, pd.Timedelta):
            is_valid = duration > pd.Timedelta(hours=2)
        else:
            is_valid = duration > 120
            
        # Spatial filter (> 2 degrees)
        lon_diff = group['longitude'].iloc[-1] - group['longitude'].iloc[0]
        lat_diff = group['latitude'].iloc[-1] - group['latitude'].iloc[0]
        spatial_dist = np.sqrt(lon_diff**2 + lat_diff**2)
        
        is_valid = is_valid and (spatial_dist > 2.0)
            
        if is_valid:
            valid_cells.append(cell_id)

    # 4. GENERATE UNIQUE COLORS & GET DYNAMIC EXTENT
    num_tracks = len(valid_cells)

    if num_tracks == 0:
        print("No valid tracks found matching the criteria!")
        ds.close()
        return

    cmap = plt.get_cmap('turbo')
    track_colors = [cmap(i) for i in np.linspace(0, 1, num_tracks)]

    df_valid = df_clean[df_clean['cell'].isin(valid_cells)]
    lon_min, lon_max = df_valid['longitude'].min(), df_valid['longitude'].max()
    lat_min, lat_max = df_valid['latitude'].min(), df_valid['latitude'].max()

    # 5. Set up the map with Cartopy
    fig = plt.figure(figsize=(12, 8))
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())

    padding = 2.0
    ax.set_extent([lon_min - padding, lon_max + padding, 
                   lat_min - padding, lat_max + padding], crs=ccrs.PlateCarree())

    ax.coastlines(linewidth=1)
    ax.add_feature(cfeature.BORDERS, linestyle='-', edgecolor='black')
    ax.add_feature(cfeature.STATES, linestyle='--', edgecolor='gray')

    gl = ax.gridlines(draw_labels=True, linestyle='--', alpha=0.5)
    gl.top_labels = False
    gl.right_labels = False

    summary_lines = [f"Feature Track Endpoints for Event: {event_name} (> 2 Hours & > 2° Movement)", 
                     "=" * 70]

    # 6. PLOT TRACKS & COLLECT ENDPOINTS
    for idx, (cell_id, color) in enumerate(zip(valid_cells, track_colors)):
        track_num = idx + 1
        group = df_valid[df_valid['cell'] == cell_id].sort_values('time')
        
        lons = group['longitude'].values
        lats = group['latitude'].values
        
        # Ensure times are explicitly marked as UTC
        times = pd.to_datetime(group['time'])
        if times.dt.tz is None:
            times = times.dt.tz_localize('UTC')
        else:
            times = times.dt.tz_convert('UTC')
            
        start_time_str = times.iloc[0].strftime('%Y-%m-%d %H:%M:%S UTC')
        end_time_str = times.iloc[-1].strftime('%Y-%m-%d %H:%M:%S UTC')
        
        summary_lines.append(f"Track #{track_num} (Cell ID: {cell_id})")
        summary_lines.append(f"  Start: Lat {lats[0]:.2f}, Lon {lons[0]:.2f} at {start_time_str}")
        summary_lines.append(f"  End:   Lat {lats[-1]:.2f}, Lon {lons[-1]:.2f} at {end_time_str}")
        summary_lines.append("")
        
        # Draw segments with ASTRAL exact day/night checks
        for i in range(len(lons) - 1):
            point_time = times.iloc[i]
            obs = Observer(latitude=lats[i], longitude=lons[i])
            
            try:
                s = sun(obs, date=point_time.date())
                sunrise = s['sunrise']
                sunset = s['sunset']
                is_night = (point_time < sunrise) or (point_time > sunset)
            except ValueError:
                month = point_time.month
                if lats[i] > 0:  
                    is_night = False if 4 <= month <= 9 else True
                else:            
                    is_night = True if 4 <= month <= 9 else False

            linestyle = '--' if is_night else '-'
            
            ax.plot(lons[i:i+2], lats[i:i+2], 
                    linewidth=1.5, color=color, linestyle=linestyle, 
                    transform=ccrs.PlateCarree())
        
        # Markers
        ax.plot(lons[0], lats[0], marker='^', color='green', markersize=8, markeredgecolor='black', transform=ccrs.PlateCarree(), zorder=5)
        ax.plot(lons[-1], lats[-1], marker='s', color='red', markersize=6, markeredgecolor='black', transform=ccrs.PlateCarree(), zorder=5)
        
        # Add sequential number
        ax.text(lons[0] + 0.3, lats[0] + 0.3, str(track_num), color='black', 
                fontsize=10, weight='bold', transform=ccrs.PlateCarree(), zorder=6)
        ax.background_img(name='NE', resolution='highest')

    plt.title(f'Event: {event_name} | {num_tracks} Tracks (> 2 Hrs, > 2° Dist)\nSolid=Day | Dashed=Night (Local Sun) | Numbered at Start', fontsize=14, pad=10)
    plt.tight_layout()
    
    plt.savefig(output_plot, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"Plot successfully saved to: {output_plot}")

    # 7. PRINT SUMMARY TO SCREEN AND SAVE TO TEXT FILE
    with open(output_text, 'w') as file:
        for line in summary_lines:
            print(line)
            file.write(line + "\n")  

    print(f"Text summary successfully saved to: {output_text}")

    ds.close()

if __name__ == "__main__":
    main()
