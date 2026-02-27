#!/usr/bin/env python

import sys
sys.path.insert(0, '/discover/nobackup/acollow/GMAOpyobs/src')

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import yaml
from pyobs.icartt import ICARTT
import os
from glob import glob
import netCDF4 as nc
import re

yamlkey_var = sys.argv[1]
flightdate = sys.argv[2]

with open('variablemap.yaml') as f:
    var = yaml.safe_load(f)
with open('plotconfig.yaml') as f:
    config = yaml.safe_load(f)

def evaluate_geosvar_expression(expr, ds):
    """
    Evaluate a simple expression like 'OC + BR' or 'SCA / EXT'
    Supports +, -, *, / operators
    ds can be a single dataset or a dictionary of datasets
    """
    # If it's just a simple variable name (no operators), return it directly
    if not any(op in expr for op in ['+', '-', '*', '/']):
        var_name = expr.strip()
        if isinstance(ds, dict):
            # For multi-file variables, extract the base variable name
            return ds[var_name].values
        else:
            return ds[var_name].values
    
    # Find all variable names in the expression
    var_names = re.findall(r'[A-Za-z_][A-Za-z0-9_]*', expr)
    
    # Create a namespace with the variables from the dataset(s)
    if isinstance(ds, dict):
        # Multiple datasets case
        namespace = {name: ds[name].values for name in var_names if name in ds}
    else:
        # Single dataset case
        namespace = {name: ds[name].values for name in var_names if name in ds}
    
    # Evaluate the expression
    return eval(expr, {"__builtins__": {}}, namespace)

# Load observational data once
ictfile = glob(config['ictpath'] + '/asiaaq-mrg60_dc8_*' + flightdate + '*ict')
m = ICARTT(ictfile)
#print([key for key in m.__dict__.keys()])
ALT, obs_ts = m.Altitude_AGL_m_DIGANGI, m.__dict__[var[yamlkey_var]['obsvar']]

# Apply unit conversion if specified in YAML
if 'obs_unit_conversion' in var[yamlkey_var]:
    obs_ts = obs_ts * var[yamlkey_var]['obs_unit_conversion']

# Define bin edges for vertical profiles
bin_edges = np.arange(0, 4001, 250)
heightbins = np.digitize(ALT, bin_edges)
height = np.arange(0.125, 4.1, step=0.25)

# Calculate observational statistics once
obs = np.array([np.nanmean(obs_ts[heightbins == h]) if np.any(heightbins == h) 
                else np.nan for h in range(1, 17)])
obs25 = np.array([np.nanpercentile(obs_ts[heightbins == h], 25) if np.any(heightbins == h) 
                  else np.nan for h in range(1, 17)])
obs75 = np.array([np.nanpercentile(obs_ts[heightbins == h], 75) if np.any(heightbins == h) 
                  else np.nan for h in range(1, 17)])

# Store model data for each experiment
model_data = {}
colors = ['red', 'blue', 'green', 'orange', 'purple', 'brown']  # Add more colors as needed

# Loop through experiments
for exp_idx, exp in enumerate(config['experiments']):
    exp_name = exp['name']
    exp_path = exp['path']
    
    # Check if sampling is a dictionary (multi-file case) or a string (single file)
    sampling_config = var[yamlkey_var]['sampling']
    
    if isinstance(sampling_config, dict):
        # Multi-file case: load multiple files and create a combined dataset dictionary
        ds_dict = {}
        for var_alias, sampling_name in sampling_config.items():
            sampledfile = glob(exp_path + '/asiaaq-mrg60_dc8_*' + flightdate + '*' + 
                              sampling_name + '.nc4')
            
            if not sampledfile:
                print(f"Warning: No sampled file found for {var_alias} in experiment {exp_name}")
                continue
            
            # Load the dataset
            ds_temp = xr.open_dataset(sampledfile[0], engine='netcdf4')
            
            # Store altitude from first file
            if 'alt_geos' not in locals():
                alt_geos = ds_temp['H'].values
            
            # Extract the actual variable name from the alias (e.g., 'SCA_rh8' -> 'SCA')
            actual_var = var_alias.split('_')[0]
            ds_dict[var_alias] = ds_temp[actual_var]
            ds_temp.close()
        
        # Evaluate expression using the dictionary of datasets
        vardata_geos = evaluate_geosvar_expression(var[yamlkey_var]['geosvar'], ds_dict)
        
    else:
        # Single file case (original behavior)
        sampledfile = glob(exp_path + '/asiaaq-mrg60_dc8_*' + flightdate + '*' + 
                          sampling_config + '.nc4')
        
        if not sampledfile:
            print(f"Warning: No sampled file found for experiment {exp_name}")
            continue
        
        # Load data
        ds = xr.open_dataset(sampledfile[0], engine='netcdf4')
        alt_geos = ds['H'].values
        
        # Evaluate the geosvar expression (handles both simple variables and expressions)
        vardata_geos = evaluate_geosvar_expression(var[yamlkey_var]['geosvar'], ds)
        
        ds.close()
    
    # Get the indices of minimum distances for all columns at once
    indices = np.array([np.argmin(np.abs(alt_geos[n, :] - ALT[n])) for n in range(len(ALT))])
    
    # Use advanced indexing to get the corresponding values from GEOS and filter for available obs
    geos_ts = np.array([vardata_geos[n, indices[n]] for n in range(len(ALT))])
    geos_ts[np.isnan(obs_ts)] = np.nan
    
    # Find layer stats for this experiment
    geos = np.array([np.nanmean(geos_ts[heightbins == h]) if np.any(heightbins == h) 
                     else np.nan for h in range(1, 17)])
    geos25 = np.array([np.nanpercentile(geos_ts[heightbins == h], 25) if np.any(heightbins == h) 
                       else np.nan for h in range(1, 17)])
    geos75 = np.array([np.nanpercentile(geos_ts[heightbins == h], 75) if np.any(heightbins == h) 
                       else np.nan for h in range(1, 17)])
    
    # Store results
    model_data[exp_name] = {
        'mean': geos,
        'p25': geos25,
        'p75': geos75,
        'color': exp.get('color', colors[exp_idx % len(colors)])
    }

# Plot
fig = plt.figure(figsize=(10, 8))
fig.tight_layout(pad=1.0)
plt.rcParams['font.size'] = '14'
ax2 = plt.subplot(111)

# Plot observations
plt.plot(obs, height, 'black', linewidth=2, label='Obs')
ax2.fill_betweenx(height, obs25, obs75, facecolor='dimgray', alpha=0.3, 
                   label='Obs 25th-75th')

# Plot each model experiment
for exp_name, data in model_data.items():
    plt.plot(data['mean'], height, color=data['color'], linewidth=2, label=exp_name)
    ax2.fill_betweenx(height, np.squeeze(data['p25']), np.squeeze(data['p75']), 
                      facecolor=data['color'], alpha=0.2, 
                      label=f"{exp_name} 25th-75th")

plt.xlabel(var[yamlkey_var]['fullname'] + ' ' + var[yamlkey_var]['units'])
plt.ylabel('Height (km)')
ax2.legend(fontsize=10, loc='best')

country = config['flightdate'][int(flightdate)]['country']
datestr = config['flightdate'][int(flightdate)]['datestring']
plt.title(f"{country}, {datestr}")

# Save with all experiment names in filename
exp_names = '_'.join([exp['name'] for exp in config['experiments']])
plt.savefig(f"{exp_names}_{yamlkey_var}_{flightdate}.png", 
            dpi=300, bbox_inches='tight')
plt.show()
