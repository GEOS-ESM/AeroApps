import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import seaborn as sns
import numpy as np
import os

model = 'cam'
input_csv = f'/discover/nobackup/caturne4/AeroApps/src/Components/modules/model_utils/post/final_analysis/{model}/{model}_aeronet_fcst_analysis_timeseries_hourly_20260529_20260721.csv'
output_dir = f'./chiclet_plots/{model}_data/'  # Change this to your desired output directory

os.makedirs(output_dir, exist_ok=True)

print(f"Loading data from {input_csv}...")
df = pd.read_csv(input_csv)

stations = ['GSFC', 'NEON_UNDE','Granite_Island','Chiwaukee_Prairie','NPU_Chicago_IL',
            'London-CDN','Egbert','Toronto', 'Dayton', 'East_Brunswick', 'Flushing','UMBC','Billerica']
df = df[df['station'].isin(stations)]


start_time = "2026-07-15 00:00:00"


start_time = pd.to_datetime(start_time)
df['time'] = pd.to_datetime(df['time'])
df['initialization time'] = pd.to_datetime(df['initialization time'])
df = df[df['initialization time'] > start_time]

stations = df['station'].unique()
print(f"Found {len(stations)} unique stations. Generating plots...")


for station in stations:

    stn_df = df[df['station'] == station].copy()
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
    combined_df = combined_df.replace([-9999, -999, -999.99, 0.0, -99], np.nan)

    

    fig, ax = plt.subplots(figsize=(12, 6))
    sns.heatmap(combined_df, 
                cmap='YlOrRd', 
                cbar_kws={'label': 'AOD'},
                linewidths=0.5, 
                linecolor='lightgray',
                ax=ax)

    times = combined_df.columns
    x_labels = [t.strftime('%m-%d %Hz') if pd.notnull(t) else '' for t in times]
    

    tick_frequency = 6

    
    ax.set_xticks(np.arange(0.5, len(x_labels), tick_frequency))
    ax.set_xticklabels([x_labels[i] for i in range(0, len(x_labels), tick_frequency)], rotation=45, ha='right')
    
    ax.set_title(f'AOD Forecast Evolution in GEOS-{model.upper()}\nStation: {station}',loc = 'left', fontsize=14, pad=15)
    ax.set_title(f'Lat: {lat:.2f}, Lon: {lon:.2f}', loc = 'right', fontsize = 14, pad = 15)
    ax.set_xlabel('Valid Time', fontsize=12)
    ax.set_ylabel('Data Source / Forecast Init', fontsize=12)

    ax.invert_yaxis()
    
    plt.tight_layout()

    out_file = os.path.join(output_dir, f'chiclet_aod_{station}.png')
    plt.savefig(out_file, dpi=300, bbox_inches='tight')
    plt.close()

print(f"Plots have been saved to {output_dir}")
