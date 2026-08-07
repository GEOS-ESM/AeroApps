import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import seaborn as sns
import numpy as np
import os
import warnings
warnings.filterwarnings('ignore')


def plot_chiclet(model, data, axs, station, lat, lon, vmin, vmax):
    # --- Plotting ---
    if model == 'cam':
        ax = 0
    elif model =='fp':
        ax = 1

    sns.heatmap(data, 
                cmap='YlOrRd', 
                vmin=vmin,
                vmax=vmax,
                cbar=False, 
                linewidths=0.5, 
                linecolor='lightgray',
                ax=axs[ax])
    
    # Format the X-axis (Time)
    times = data.columns
    x_labels = [t.strftime('%m-%d %Hz') if pd.notnull(t) else '' for t in times]
    
    # Show every 6th tick
    tick_frequency = 6 
    axs[ax].set_xticks(np.arange(0.5, len(x_labels), tick_frequency))
    axs[ax].set_xticklabels([x_labels[i] for i in range(0, len(x_labels), tick_frequency)], rotation=45, ha='right')
    axs[ax].set_title(f'AOD Forecast Evolution in GEOS-{model.upper()}\nStation: {station}', loc='left', fontsize=20, pad=15)
    axs[ax].set_title(f'Lat: {lat:.2f}, Lon: {lon:.2f}', loc='right', fontsize=20, pad=15)
    axs[ax].set_xlabel('Valid Time', fontsize=18)
    axs[ax].set_ylabel('Data Source / Forecast Init', fontsize=18)
    axs[ax].tick_params(axis ='both', labelsize = 16)
    axs[ax].invert_yaxis()


# --- Configuration ---
models = ['cam', 'fp']
stations = ['GSFC']

output_dir = './chiclet_plots/combined_data/'
os.makedirs(output_dir, exist_ok=True)

model_data = {}
for model in models:
    input_csv = f'/discover/nobackup/caturne4/AeroApps/src/Components/modules/model_utils/post/final_analysis/{model}/{model}_aeronet_fcst_analysis_timeseries_hourly_20260529_20260721.csv'
    
    print(f"Loading data from {input_csv}...")
    df = pd.read_csv(input_csv)
    df = df[df['station'].isin(stations)]

    custom_period = True
    if custom_period == True:
        start_time = "2026-07-15 12:00:00"
        end_time   = "2026-07-20 00:00:00"
        start_time = pd.to_datetime(start_time)
        end_time   = pd.to_datetime(end_time)
        
        df['time'] = pd.to_datetime(df['time'])
        df['initialization time'] = pd.to_datetime(df['initialization time'])
        df = df[df['initialization time'] >= start_time]
        df = df[df['initialization time'] <= end_time]
        df = df[df['time'] <= end_time]
        output_dir = './chiclet_plots/custom_period/'

    
    model_data[model] = df
    
unique_stations = model_data['cam']['station'].unique()
print(f"Found {len(unique_stations)} unique stations. Generating plots...")

# Process each station
for station in unique_stations:
    station_plot_data = {}
    global_min = np.inf
    global_max = -np.inf
    
    for model in models:
        current_model_df = model_data[model]
        stn_df = current_model_df[current_model_df['station'] == station].copy()
        
        if stn_df.empty:
            continue
            
        stn_df = stn_df.sort_values(by='time')
        lat = stn_df['lats'].min()
        lon = stn_df['lons'].min()
        obs_df = stn_df.drop_duplicates('time')[['time', 'aeronet_aod']].set_index('time')
        obs_df.columns = ['AERONET Obs']
        ana_df = stn_df.drop_duplicates('time')[['time', 'analysis_aod']].set_index('time')
        ana_df.columns = [f'GEOS-{model.upper()} Analysis']
        
        init_times = sorted(stn_df['initialization time'].unique())
        fcst_dict = {}
        for init in init_times:
            init_df = stn_df[stn_df['initialization time'] == init][['time', 'model_aod']].set_index('time')
            row_label = f"Init: {pd.to_datetime(init).strftime('%m/%d %Hz')}"
            fcst_dict[row_label] = init_df['model_aod']
        
        fcst_df = pd.DataFrame(fcst_dict)
        combined_df = pd.concat([obs_df, ana_df, fcst_df], axis=1).T
        station_plot_data[model] = {'data': combined_df, 'lat': lat, 'lon': lon}
        
        
        global_min = min(global_min, combined_df.min().min())
        global_max = max(global_max, combined_df.max().max())

    if not station_plot_data:
        continue
   
    fig, axs = plt.subplots(2, 1, figsize=(18, 12))
    for model in models:
        plot_data = station_plot_data[model]
        plot_chiclet(model, plot_data['data'], axs, station, plot_data['lat'], plot_data['lon'], global_min, global_max)
        
 
    sm = plt.cm.ScalarMappable(cmap='YlOrRd', norm=plt.Normalize(vmin=global_min, vmax=global_max))
    sm.set_array([]) # Required for matplotlib 3.1 and later
    cbar_ax = fig.add_axes([1.02, 0.15, 0.02, 0.7]) 
    cbar = fig.colorbar(sm, cbar_ax, pad = 0.1)
    cbar.set_label('AOD', fontsize=18)
    cbar.ax.tick_params(labelsize = 16)
    fig.tight_layout()
    

    out_file = os.path.join(output_dir, f'chiclet_aod_{station}.png')
    plt.savefig(out_file, dpi=300, bbox_inches='tight')
    plt.close()
    
print(f"Done! Plots have been saved to {output_dir}")