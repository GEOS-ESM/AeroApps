'''
Plots field campaign profiles of organic carbon, sulfate, and nitrate from already
sampled files from MERRA-21C and MERRA-2 in comparison to AMS observations.

export SRC_DIR="/gpfsm/dnb34/acollow/AeroApps/AeroApps"
export PYTHONPATH="${SRC_DIR}/install/lib/Python"
'''

import os
import glob
import yaml
import string
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from pyobs.icartt import ICARTT

plt.rcParams.update({
    'font.size': 20,
    'axes.labelsize': 22,
    'axes.titlesize': 24,
    'xtick.labelsize': 20,
    'ytick.labelsize': 20,
    'legend.fontsize': 18,
    'lines.linewidth': 2.5
})

# Define the campaigns and their corresponding YAML config files
campaigns = [
    {'name': 'SEAC4RS', 'yaml': 'sampling_seac4rsBB.yaml'},
    {'name': 'FIREX-AQ', 'yaml': 'sampling_firex.yaml'},
    {'name': 'CAMP2Ex', 'yaml': 'sampling_camp2ex.yaml'},
    {'name': 'ASIA-AQ (Thailand)', 'yaml': 'sampling_asiaaqthailand.yaml'}
]

species_list = ['oc', 'sulfate', 'nitrate']
species_titles = ['Organic Carbon', 'Sulfate', 'Nitrate']

ALT_BINS = np.arange(0, 10.5, 0.5)
ALT_CENTERS = (ALT_BINS[:-1] + ALT_BINS[1:]) / 2.0

def get_binned_stats(altitude, values):
    """Calculate median, 25th, and 75th percentiles for altitude bins."""
    medians = np.full(len(ALT_CENTERS), np.nan)
    p25 = np.full(len(ALT_CENTERS), np.nan)
    p75 = np.full(len(ALT_CENTERS), np.nan)
    
    # Filter out NaNs
    valid = ~np.isnan(altitude) & ~np.isnan(values)
    alt_valid = altitude[valid]
    val_valid = values[valid]
    
    for i in range(len(ALT_BINS)-1):
        idx = (alt_valid >= ALT_BINS[i]) & (alt_valid < ALT_BINS[i+1])
        if np.sum(idx) > 5:
            bin_data = val_valid[idx]
            medians[i] = np.nanmedian(bin_data)
            p25[i] = np.nanpercentile(bin_data, 25)
            p75[i] = np.nanpercentile(bin_data, 75)
            
    return medians, p25, p75

def extract_campaign_data(config, species_key):
    """Extracts and concatenates arrays for obs, M21C, and M2, masking models to valid obs."""

    obs_vals, m21c_vals, m2_vals = [], [], []
    obs_alts, m21c_alts, m2_alts = [], [], []

    obs_var = config[species_key]['obs']
    m21c_var = config[species_key]['m21c']
    
    # MERRA-2 does not have nitrate
    m2_var = config.get(species_key, {}).get('m2', None)
    if species_key == 'nitrate':
        m2_var = None 

    alt_var_obs = config['altitude']['obs']
    alt_units_obs = str(config['altitude'].get('obsunits', 'm')).strip().lower()
    alt_var_m21c = config['altitude']['m21c']

    merge_dir = config['mergefiles']
    out_dir = config['sampled_outdir']

    for date in config.get('flight_dates', []):
        print(date)
        obs_step_vals = None
        obs_step_alts = None
        valid_time_mask = None
        current_m21c_alt = None
        matched_m21c_vals = None

        # ---------------------------------------------------------
        # 1. Observations (ICARTT)
        # ---------------------------------------------------------
        obs_files = glob.glob(os.path.join(merge_dir, f"*{date}*.ict"))
        if obs_files:
            try:
                obs_data = ICARTT(obs_files[0], to_NaN=True)
                
                if hasattr(obs_data, obs_var) and hasattr(obs_data, alt_var_obs):
                    obs_step_vals = getattr(obs_data, obs_var)
                    obs_step_alts = getattr(obs_data, alt_var_obs)
                elif type(obs_data) == dict and obs_var in obs_data and alt_var_obs in obs_data:
                    obs_step_vals = obs_data[obs_var]
                    obs_step_alts = obs_data[alt_var_obs]
                elif hasattr(obs_data, 'df') and obs_var in obs_data.df.columns and alt_var_obs in obs_data.df.columns:
                    obs_step_vals = obs_data.df[obs_var].values
                    obs_step_alts = obs_data.df[alt_var_obs].values
                
                if obs_step_vals is not None:
                    valid_time_mask = ~np.isnan(obs_step_vals) & ~np.isnan(obs_step_alts)
                    
                    obs_vals.extend(obs_step_vals[valid_time_mask])
                    obs_alts.extend(obs_step_alts[valid_time_mask])
            except Exception:
                pass

        # ---------------------------------------------------------
        # 2. MERRA-21C
        # ---------------------------------------------------------
        if out_dir.endswith('MERRA-21C'):
            m21c_search = os.path.join(out_dir, f"*{date}*.pm_ams.nc4")
        else:
            m21c_search = os.path.join(out_dir, 'MERRA-21C', f"*{date}*.pm_ams.nc4")
            
        m21c_files = glob.glob(m21c_search)
        if m21c_files and valid_time_mask is not None:
            try:
                ds = xr.open_dataset(m21c_files[0])
                if '+' in m21c_var:
                    v_raw = sum([ds[p].values for p in m21c_var.split('+')])
                else:
                    v_raw = ds[m21c_var].values
                a_raw = ds[alt_var_m21c].values
                ds.close()
                
                if len(v_raw) == len(valid_time_mask):
                    v_raw_valid = v_raw[valid_time_mask]
                    a_raw_valid = a_raw[valid_time_mask]
                    obs_alts_valid = obs_step_alts[valid_time_mask]

                    if v_raw_valid.ndim > 1:
                        abs_diff = np.abs(a_raw_valid - obs_alts_valid[:, None])
                        closest_level_idx = np.argmin(abs_diff, axis=1)
                        
                        matched_m21c_vals = v_raw_valid[np.arange(len(v_raw_valid)), closest_level_idx]
                        matched_m21c_alts = a_raw_valid[np.arange(len(a_raw_valid)), closest_level_idx]
                    else:
                        matched_m21c_vals = v_raw_valid
                        matched_m21c_alts = a_raw_valid

                    m21c_vals.extend(matched_m21c_vals)
                    m21c_alts.extend(matched_m21c_alts)
                    current_m21c_alt = closest_level_idx if v_raw_valid.ndim > 1 else None
            except Exception as e:
                print(f"  -> WARNING: Error processing MERRA-21C for {date}: {e}")
        elif not m21c_files:
            print(f"  -> WARNING: Missing MERRA-21C file for date {date}")

        # ---------------------------------------------------------
        # 3. MERRA-2
        # ---------------------------------------------------------
        if m2_var is not None and valid_time_mask is not None:
            if 'MERRA-21C' in out_dir:
                m2_dir = out_dir.replace('MERRA-21C', 'MERRA-2')
            else:
                m2_dir = os.path.join(out_dir, 'MERRA-2')
                
            m2_search = os.path.join(m2_dir, f"*{date}*.pm_ams.nc4")
            m2_files = glob.glob(m2_search)
            
            if m2_files:
                try:
                    ds = xr.open_dataset(m2_files[0])
                    v_raw = ds[m2_var].values
                    ds.close()
                    
                    if len(v_raw) == len(valid_time_mask):
                        v_raw_valid = v_raw[valid_time_mask]
                        
                        if v_raw_valid.ndim > 1 and current_m21c_alt is not None:
                            matched_m2_vals = v_raw_valid[np.arange(len(v_raw_valid)), current_m21c_alt]
                        else:
                            matched_m2_vals = v_raw_valid

                        m2_vals.extend(matched_m2_vals)
                        m2_alts.extend(matched_m21c_alts if 'matched_m21c_alts' in locals() else [])
                except Exception:
                    pass
    
    # ---------------------------------------------------------
    # Put all altitudes in km
    # ---------------------------------------------------------
    obs_alts_arr = np.array(obs_alts)
    obs_vals_arr = np.array(obs_vals)
    
    m21c_alts_arr = np.array(m21c_alts) / 1000.0 if len(m21c_alts) > 0 else np.array([])
    m21c_vals_arr = np.array(m21c_vals)
    
    m2_alts_arr = np.array(m2_alts) / 1000.0 if len(m2_alts) > 0 else np.array([])
    m2_vals_arr = np.array(m2_vals)
    
    if alt_units_obs == 'm' and len(obs_alts_arr) > 0:
        obs_alts_arr = obs_alts_arr / 1000.0

    return {
        'obs': (obs_alts_arr, obs_vals_arr),
        'm21c': (m21c_alts_arr, m21c_vals_arr),
        'm2': (m2_alts_arr, m2_vals_arr) if m2_var else None
    }

# -----------------
# Generate the Plot
# -----------------
print("Generating profiles... Please wait.")
fig, axes = plt.subplots(nrows=4, ncols=3, figsize=(18, 22), sharey=True)
panel_labels = list(string.ascii_lowercase)
panel_idx = 0

for row, campaign_info in enumerate(campaigns):
    with open(campaign_info['yaml'], 'r') as f:
        config = yaml.safe_load(f)
    
    for col, (species, title) in enumerate(zip(species_list, species_titles)):
        ax = axes[row, col]
        
        data = extract_campaign_data(config, species)
        
        colors = {'obs': 'blue', 'm21c': 'black', 'm2': 'red'}
        labels = {'obs': 'Obs', 'm21c': 'MERRA-21C', 'm2': 'MERRA-2'}
        
        keys_to_plot = ['obs', 'm21c']
        if data['m2'] is not None:
            keys_to_plot.append('m2')
            
        for key in keys_to_plot:
            alt_vals, val_vals = data[key]
            
            # Skip if arrays are empty
            if len(alt_vals) == 0:
                continue
                
            med, p25, p75 = get_binned_stats(alt_vals, val_vals)
            
            ax.plot(med, ALT_CENTERS, color=colors[key], linestyle='-', linewidth=2.5, label=labels[key])
            ax.fill_betweenx(ALT_CENTERS, p25, p75, color=colors[key], alpha=0.2, edgecolor='none')
            

        ax.set_ylim(0, 8) # 0 to 10 km
        ax.grid(True, linestyle='--', alpha=0.5)
        
        ax.text(0.005, 0.05, f'({panel_labels[panel_idx]})', transform=ax.transAxes, 
                fontsize=16, fontweight='bold', va='top', ha='left')
        panel_idx += 1
        
        if col == 0:
            ax.set_ylabel(f"{campaign_info['name']}\nAltitude (km)")
        if row == 0:
            ax.set_title(title, fontweight='bold')
        if row == 3:
            ax.set_xlabel(r'Mass Concentration ($\mu g/m^3$)')
            #ax.set_ylim(0,4)
        if row == 2:
            ax.set_ylim(0,4)
        if row == 0:
            if col == 0:
                ax.legend(loc='upper right')

plt.tight_layout()
plt.subplots_adjust(top=0.88, hspace=0.15, wspace=0.1)

# Save the figure
plt.savefig('smoke_AMSprofiles_campaigns.png', dpi=300, bbox_inches='tight')
print("Plot successfully generated and saved as smoke_AMSprofiles_campaigns.png")
plt.show()
