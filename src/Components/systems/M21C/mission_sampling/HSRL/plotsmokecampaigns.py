import os
import yaml
import glob
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from scipy.interpolate import griddata, interp1d
import warnings

warnings.filterwarnings('ignore')

# Set publication quality plot parameters
plt.rcParams.update({
    'font.size': 16,
    'axes.labelsize': 18,
    'axes.titlesize': 20,
    'xtick.labelsize': 16,
    'ytick.labelsize': 16,
    'legend.fontsize': 16,
    'lines.linewidth': 2.5
})

def load_config(yaml_file):
    with open(yaml_file, 'r') as f:
        return yaml.safe_load(f)

def get_obs_var_name(model_var, wavelength, instrument):
    """Maps the model variable to the exact observation variable name format based on instrument."""
    instrument = str(instrument).upper()
    
    if instrument == "DIAL":
        mapping = {
            "EXT": f"ext_{wavelength}nm_prfl",
            "SCA": f"bsc_{wavelength}nm_prfl", 
            "BSC": f"bsc_{wavelength}nm_prfl", 
            "DEPOL": f"aerdep_{wavelength}nm_prfl"
        }
    else:
        # HALO and HSRL2
        mapping = {
            "EXT": f"{wavelength}_ext", 
            "SCA": f"{wavelength}_bsc_cloud_screened", 
            "BSC": f"{wavelength}_bsc_cloud_screened", 
            "DEPOL": f"{wavelength}_aer_dep"
        }
        
    return mapping.get(model_var, f"{wavelength}_{model_var.lower()}")

def standardize_time_to_hours(time_array, flight_date_str):
    """
    Standardizes any time array (datetime64, decimal hours, or seconds) 
    into Decimal Hours since midnight UTC of the flight date.
    """
    date_str = str(flight_date_str)
    # Define midnight UTC of the flight date
    epoch = np.datetime64(f"{date_str[:4]}-{date_str[4:6]}-{date_str[6:8]}T00:00:00")
    
    if np.issubdtype(time_array.dtype, np.datetime64):
        # Convert datetime64 directly to hours since epoch
        return (time_array - epoch) / np.timedelta64(1, 's') / 3600.0
    else:
        time_float = time_array.astype(float)
        if np.nanmax(time_float) < 48.0:
            # It's already in decimal hours (e.g., typical gps_time)
            return time_float
        elif np.nanmax(time_float) > 1e8:
            # It's in Unix Epoch seconds
            epoch_unix = epoch.astype('datetime64[s]').astype(float)
            return (time_float - epoch_unix) / 3600.0
        else:
            # It's in seconds since midnight
            return time_float / 3600.0

def process_flight_data(config, date, target_variable, model_versions, wavelength):
    """
    Loads and processes data for a single flight for multiple models. 
    Handles multiple flight legs (e.g., _L1, _L2) per day seamlessly.
    """
    model_dir = config.get('output_dir', '')
    obs_dir = config.get('obs_dir', '')
    instrument = config.get('instrument', 'HALO').upper()
    
    # Try different cases for the campaign name
    campaign_original = config['campaign']
    campaign_variations = [
        campaign_original,
        campaign_original.upper(),
        campaign_original.lower()
    ]
    
    obs_files = []
    for camp_name in campaign_variations:
        search_pattern = f"{camp_name}-{instrument}_{config['plane']}_{date}*.h5"
        matches = glob.glob(os.path.join(obs_dir, search_pattern))
        if matches:
            obs_files.extend(matches)
            
    obs_files = sorted(list(set(obs_files)))
            
    if not obs_files:
        print(f"Missing OBS data for date {date}. Skipping...")
        return None
        
    all_points = []
    all_values = []
    b_kwargs = {'phony_dims': 'sort'}
    groups_to_check = [None, 'DataProducts', 'Data_Products', 'Nav_Data', 'State', 'Navigation', 'ApplanixIMU', 'UserInput']
    
    for obs_file in obs_files:
        ds_obs = xr.open_dataset(obs_file, engine='h5netcdf', backend_kwargs=b_kwargs)
        obs_var_name = get_obs_var_name(target_variable, wavelength, instrument)
        
        if obs_var_name not in ds_obs:
            found = False
            for group in ['Data_Products', 'DataProducts']:
                try:
                    ds_obs_group = xr.open_dataset(obs_file, group=group, engine='h5netcdf', backend_kwargs=b_kwargs)
                    if obs_var_name in ds_obs_group:
                        ds_obs = ds_obs_group
                        found = True
                        break
                except Exception:
                    pass
            if not found:
                continue

        obs_data = ds_obs[obs_var_name].values  
        
        # Apply DIAL specific Cloud Mask
        if instrument == "DIAL":
            cloud_mask_name = "Cloud_Mask_prfl"
            if cloud_mask_name in ds_obs:
                cloud_mask = ds_obs[cloud_mask_name].values
                if cloud_mask.shape == obs_data.shape:
                    obs_data = np.where(cloud_mask == 1, np.nan, obs_data)

        # Mask out negative noise, NaNs, and unphysically large fill values (e.g. 999.0)
        obs_data = np.where((obs_data < 0) | (obs_data > 100) | np.isnan(obs_data), np.nan, obs_data)
        
        # --- NEW: Remove points > 2 standard deviations from the mean ---
        obs_mean = np.nanmean(obs_data)
        obs_std = np.nanstd(obs_data)
        obs_data = np.where(np.abs(obs_data - obs_mean) > (2 * obs_std), np.nan, obs_data)
        # ----------------------------------------------------------------

        obs_time_var = None
        for group in groups_to_check:
            try:
                ds_temp = xr.open_dataset(obs_file, group=group, engine='h5netcdf', backend_kwargs=b_kwargs)
                for k in ['time', 'Midtime', 'gps_time', 'Time', 'Time_UTC']:
                    if k in ds_temp:
                        obs_time_var = ds_temp[k]
                        break
                if obs_time_var is not None:
                    break
            except Exception:
                pass

        obs_alt_var = None
        for group in groups_to_check:
            try:
                ds_temp = xr.open_dataset(obs_file, group=group, engine='h5netcdf', backend_kwargs=b_kwargs)
                # gps_alt removed so it doesn't grab the 1D aircraft altitude array
                for k in ['Altitude', 'Altitudes', 'z']:
                    if k in ds_temp:
                        obs_alt_var = ds_temp[k]
                        break
                if obs_alt_var is not None:
                    break
            except Exception:
                pass

        if obs_time_var is None or obs_alt_var is None:
            continue

        obs_time = obs_time_var.values.squeeze()
        obs_lev = obs_alt_var.values.squeeze()
        
        # Explicitly convert Altitudes from meters to kilometers
        obs_lev = obs_lev / 1000.0  

        obs_time_hours = standardize_time_to_hours(obs_time, date)

        if obs_data.shape == (len(obs_lev), len(obs_time)):
            aligned_obs_data = obs_data          
        elif obs_data.shape == (len(obs_time), len(obs_lev)):
            aligned_obs_data = obs_data.T        
        else:
            aligned_obs_data = obs_data          

        if obs_time_hours.ndim == 1 and obs_lev.ndim == 1:
            T_obs, Z_obs = np.meshgrid(obs_time_hours, obs_lev)
        else:
            T_obs, Z_obs = obs_time_hours, obs_lev

        # Flatten points for interpolation
        T_flat = T_obs.flatten()
        Z_flat = Z_obs.flatten()
        V_flat = aligned_obs_data.flatten()

        # 4. Check for Surface Elevation to convert MSL to AGL
        surface_elev_var = None
        sfc_varnames = [
            'DEM_gnd_alt', 'DEM_elevation', 'Surface_Elevation', 
            'Terrain_Elevation', 'Surface_Alt', 'DEM_altitude', 
            'DEM_alt', 'DEM_Alt', 'elevation', 'Elevation', 'surface_alt'
        ]
        
        for group in groups_to_check:
            try:
                ds_temp = xr.open_dataset(obs_file, group=group, engine='h5netcdf', backend_kwargs=b_kwargs)
                for k in sfc_varnames:
                    if k in ds_temp:
                        surface_elev_var = ds_temp[k]
                        break
                if surface_elev_var is not None:
                    break
            except Exception:
                pass

        if surface_elev_var is not None:
            # Squeeze to 1D time array
            sfc_elev = surface_elev_var.values.squeeze()
            
            # Mask out missing data fill values (e.g. -999.0 or -9999.0) so they don't skew the AGL math
            sfc_elev = np.where(sfc_elev < -900, np.nan, sfc_elev)
            
            # Explicitly convert DEM_altitude from meters to kilometers
            sfc_elev = sfc_elev / 1000.0
                
            # Use interpolation to map the 1D surface elevation exactly to the flattened Time points
            valid_sfc = ~np.isnan(sfc_elev)
            if np.any(valid_sfc):
                sfc_interp_func = interp1d(obs_time_hours[valid_sfc], sfc_elev[valid_sfc], 
                                           bounds_error=False, fill_value="extrapolate")
                
                # Subtract surface elevation from altitude points to create AGL
                sfc_at_points = sfc_interp_func(T_flat)
                Z_flat = Z_flat - sfc_at_points
                
                # Mask out underground data
                V_flat = np.where(Z_flat < 0, np.nan, V_flat)
        else:
            pass # No terrain found, will remain in ASL.

        # Filter out NaNs before appending
        valid_mask = ~np.isnan(V_flat)
        
        leg_points = np.array([T_flat[valid_mask], Z_flat[valid_mask]]).T
        leg_values = V_flat[valid_mask]
        
        all_points.append(leg_points)
        all_values.append(leg_values)
        
    if not all_points or sum(len(v) for v in all_values) == 0:
        return None

    points = np.vstack(all_points)
    values = np.concatenate(all_values)

    flight_results = {'models': {}, 'obs': None, 'lev': None}
    collection = f"ext{wavelength}"
    
    aer_file_pattern = f"{config['campaign']}-MERRA21C-aer-Nv-{config['plane']}_Model_{date}_R*.nc*"
    aer_files = glob.glob(os.path.join(model_dir, aer_file_pattern))
    
    if not aer_files:
        return None
        
    ds_aer = xr.open_dataset(aer_files[0])
    if 'H' not in ds_aer:
        return None
        
    model_H = ds_aer['H'].values
    if np.nanmax(model_H) > 100:  
        Z_mod_km = model_H / 1000.0
    else:
        Z_mod_km = model_H

    # --- NEW FIX: Convert Model H from ASL to AGL ---
    # Try to find the model's surface height (ZS)
    if 'ZS' in ds_aer:
        zs = ds_aer['ZS'].values
        if np.nanmax(zs) > 100:
            zs = zs / 1000.0
        # Subtract surface elevation
        if zs.ndim == 1:
            Z_mod_km = Z_mod_km - zs[:, None]
        elif zs.ndim == 2:
            Z_mod_km = Z_mod_km - zs
    else:
        # Fallback if ZS isn't available: subtract the lowest level 
        # of H at each time step to approximate 0 km AGL
        z_min = np.nanmin(Z_mod_km, axis=1, keepdims=True)
        Z_mod_km = Z_mod_km - z_min
        
    # Ensure no negative values due to slight interpolation anomalies
    Z_mod_km = np.where(Z_mod_km < 0, 0, Z_mod_km)
    # ------------------------------------------------
        
    # Sort model heights ascending to prevent flipped interpolations
    sort_idx = np.argsort(np.nanmean(Z_mod_km, axis=0))
    Z_mod_km = Z_mod_km[:, sort_idx]

    for mv in model_versions:
        model_file_pattern = f"{config['campaign']}-{mv}-{collection}-{config['plane']}_Model_{date}_R*.nc*"
        model_files = glob.glob(os.path.join(model_dir, model_file_pattern))
        
        if not model_files:
            continue
            
        ds_model = xr.open_dataset(model_files[0])
        
        if target_variable not in ds_model:
            continue
            
        model_data = ds_model[target_variable].values
        
        # Apply the same sorting to the model data array to match the sorted Z_mod_km
        if model_data.shape[1] == len(sort_idx):
            model_data = model_data[:, sort_idx]
            
        model_time = ds_model['time'].values

        model_time_hours = standardize_time_to_hours(model_time, date)
        
        T_mod = np.tile(model_time_hours[:, None], (1, Z_mod_km.shape[1]))
        obs_interp = griddata(points, values, (T_mod, Z_mod_km), method='linear')
        model_masked = np.where(np.isnan(obs_interp), np.nan, model_data)
        
        flight_results['models'][mv] = model_masked
        if flight_results['obs'] is None:
            flight_results['obs'] = obs_interp
            flight_results['lev'] = np.nanmean(Z_mod_km, axis=0) 
            
    return flight_results

def get_campaign_profile_data(config_path, target_variable="EXT", model_versions=["MERRA2", "MERRA21C"], wavelength=532):
    """Extracts aggregate stats for an entire campaign."""
    config = load_config(config_path)
    flight_dates = config.get('flight_dates', [])
    
    all_models_data = {mv: [] for mv in model_versions}
    all_obs_data = []
    common_lev = None
    
    for date in flight_dates:
        results = process_flight_data(config, date, target_variable, model_versions, wavelength)
        if results is not None:
            common_lev = results['lev']
            all_obs_data.append(results['obs'])
            for mv in model_versions:
                if mv in results['models']:
                    all_models_data[mv].append(results['models'][mv])
            
    if not all_obs_data:
        print(f"No valid data processed for {config_path}.")
        return None

    obs_2d = np.vstack(all_obs_data)
    obs_median = np.nanmedian(obs_2d, axis=0)
    obs_p15, obs_p85 = np.nanpercentile(obs_2d, [15, 85], axis=0)
    
    model_stats = {}
    for mv in model_versions:
        if all_models_data[mv]:
            mod_2d = np.vstack(all_models_data[mv])
            model_stats[mv] = {
                'median': np.nanmedian(mod_2d, axis=0), 
                'p15': np.nanpercentile(mod_2d, 15, axis=0),
                'p85': np.nanpercentile(mod_2d, 85, axis=0)
            }
            
    return {
        'obs_median': obs_median, 'obs_p15': obs_p15, 'obs_p85': obs_p85,
        'model_stats': model_stats, 'common_lev': common_lev,
        'instrument': config.get("instrument", "LIDAR")
    }

def main():
    target_variable = "EXT"
    wavelength = 532
    model_versions = ["MERRA2", "MERRA21C"]
    
    campaigns = [
        {"config": "config_seac4rsBB.yaml", "label": "a) SEAC4RS"},
        {"config": "config_firex.yaml", "label": "b) FIREX-AQ"},
        {"config": "config_camp2ex.yaml", "label": "c) CAMP2EX"},
        {"config": "config_oracles2018.yaml", "label": "d) ORACLES-3"}
    ]
    
    var_dict = {
        "EXT": "Extinction",
        "SCA": "Scattering",
        "BSC": "Backscatter",
        "DEPOL": "Depolarization"
    }
    var_display_name = var_dict.get(target_variable, target_variable)
    
    # Adjust titles and set conversion factor for BSC
    if target_variable == "BSC":
        x_axis_title = f"{wavelength} nm Aerosol {var_display_name} (Mm$^{{-1}}$ sr$^{{-1}}$)"
        conversion_factor = 0.001
    else:
        x_axis_title = f"{wavelength} nm Aerosol {var_display_name} (Mm$^{{-1}}$)"
        conversion_factor = 0.001
    
    colors = {"MERRA2": "red", "MERRA21C": "black"}
    
    fig, axes = plt.subplots(2, 2, figsize=(16, 16), sharey=True)
    axes = axes.flatten()
    
    for i, camp_info in enumerate(campaigns):
        print(f"Gathering data for {camp_info['label']}...")
        ax = axes[i]
        
        try:
            data = get_campaign_profile_data(camp_info['config'], target_variable, model_versions, wavelength)
        except Exception as e:
            print(f"Error processing {camp_info['config']}: {e}")
            data = None
            
        if data is None:
            ax.text(0.5, 0.5, "Data Unavailable", ha='center', va='center', transform=ax.transAxes)
            ax.set_title(camp_info['label'], loc='left', fontweight='bold')
            continue
            
        common_lev = data['common_lev']
        
        # Plot models with unit conversion
        for mv in model_versions:
            if mv in data['model_stats']:
                c = colors.get(mv, "green")
                display_name = mv.replace("MERRA", "MERRA-") if mv.startswith("MERRA") else mv
                
                mod_median = data['model_stats'][mv]['median'] / conversion_factor
                mod_p15 = data['model_stats'][mv]['p15'] / conversion_factor
                mod_p85 = data['model_stats'][mv]['p85'] / conversion_factor
                
                ax.plot(mod_median, common_lev, label=f'{display_name} Median', color=c)
                ax.fill_betweenx(common_lev, mod_p15, mod_p85, color=c, alpha=0.15)
        
        # Plot obs with unit conversion
        obs_median = data['obs_median'] / conversion_factor
        obs_p15 = data['obs_p15'] / conversion_factor
        obs_p85 = data['obs_p85'] / conversion_factor
        
        ax.plot(obs_median, common_lev, label=f'Obs Median', color='blue', linestyle='-')
        ax.fill_betweenx(common_lev, obs_p15, obs_p85, color='blue', alpha=0.15)
        
        # Axis limits formatting scaled by the conversion factor
        ax.set_ylim(0, 5)
            
        # Format the x-axis to handle scientific notation gracefully for small backscatter values
        if target_variable == "BSC":
            ax.ticklabel_format(axis='x', style='sci', scilimits=(0,0))
            
        ax.set_title(camp_info['label'], loc='left', fontweight='bold')
        ax.grid(True, linestyle=':', alpha=0.7)
        
        if i in [2, 3]: # Bottom row
            ax.set_xlabel(x_axis_title)
        if i in [0, 2]: # Left column
            ax.set_ylabel("Altitude (km)")
            
        # Only add legend to the first panel to save space and reduce clutter
        if i == 0:
            ax.legend(loc='upper right')

    plt.tight_layout()
    output_fig = f"Smoke_Campaign_Profiles_{target_variable}_{wavelength}.png"
    plt.savefig(output_fig, dpi=300, bbox_inches='tight')
    print(f"Figure saved as {output_fig}")
    plt.show()

if __name__ == "__main__":
    main()
