#!/usr/bin/env python3
"""
    Samples MERRA-2 and MERRA-21C at 60 second intervals according to flight tracks in HSRL2/DIAL/HALO files.
    Supports single or multiple legs per day.
"""
import sys
import os
import yaml
import numpy as np
from glob import glob
from datetime import datetime, timedelta
import pandas as pd
import traceback

from pyobs.sampler import TRAJECTORY
from pyobs.hsrl import HSRL
from pyobs.dial import DIAL
from dask.distributed import Client, LocalCluster

if __name__ == '__main__':
    if len(sys.argv) < 3:
        print("Usage: python3 sample_campaign.py <config.yaml> <YYYYMMDD>")
        sys.exit(1)

    config_file = sys.argv[1]
    date_str = sys.argv[2]

    with open(config_file, 'r') as f:
        config = yaml.safe_load(f)

    campaign = config['campaign']
    instrument = config['instrument']
    plane = config['plane']

    kwargs = dict(n_workers=6, threads_per_worker=1, memory_limit='4GB')
    with LocalCluster(**kwargs) as cluster, Client(cluster) as client:
        
        # Use a wildcard after the date to capture multiple legs (e.g., _L1_, _L2_, _F1_) or revisions
        h5_file_pattern_lower = f"{config['obs_dir']}/{campaign.lower()}-{instrument}_{plane}_{date_str}*.h5"
        h5_file_pattern_upper = f"{config['obs_dir']}/{campaign.upper()}-{instrument}_{plane}_{date_str}*.h5"
        
        h5_files = sorted(list(set(glob(h5_file_pattern_lower) + glob(h5_file_pattern_upper))))
        
        if not h5_files:
            print(f"Observation file missing for date {date_str}. Skipping sampling...")
            sys.exit(0)
            
        print(f"+++++++ Sampling on {len(h5_files)} file(s) for {date_str}:")
        for f in h5_files:
            print(f"  - {os.path.basename(f)}")
        
        # Determine the bracketed MERRA files ONCE for the given date
        base_dt = datetime.strptime(date_str, "%Y%m%d")
        dates_to_load = [base_dt, base_dt + timedelta(days=1)]
        
        m2_files = []
        m21c_files = []
        
        for dt_obj in dates_to_load:
            y, m, d = dt_obj.strftime("%Y"), dt_obj.strftime("%m"), dt_obj.strftime("%d")
            
            m2_pat = f"{config['merra2_archive']}/Y{y}/M{m}/MERRA2.inst3_3d_aer_Nv.{y}{m}{d}.nc4"
            m2_files.extend(glob(m2_pat))
            
            m21c_pat = f"{config['tmp_dir']}/regridded_*.{y}-{m}-{d}T*.nc4"
            m21c_files.extend(glob(m21c_pat))
            
        m2_files = sorted(m2_files)
        m21c_files = sorted(m21c_files)
        chunks = {'time': 1, 'lev': -1, 'lat': -1, 'lon': -1}

        # Loop through each flight leg/file and process them individually
        for h5_file in h5_files:
            basename = os.path.basename(h5_file)
            print(f"\n--- Processing {basename} ---")
            
            try:
                if instrument.upper() == 'DIAL':
                    if campaign.upper() == 'FIREXAQ':
                        SDS_DIAL = {
                            'Nav_Data': ('gps_alt', 'gps_lat', 'gps_lon', 'Midtime', 'DEM_gnd_alt', 'Altitudes', 'pres_alt'),
                            'DataProducts': ('O3_prfl', 'Sa_532nm_prfl', 'aerdep_532nm_prfl', 'bsc_532nm_prfl', 'depol_532nm_prfl', 'ext_532nm_prfl')
                        }
                        obs_obj = DIAL(h5_file, verbose=False, SDS=SDS_DIAL)
                    else:
                        obs_obj = DIAL(h5_file, verbose=False)
                        
                elif instrument.upper() == 'HALO':
                    from pyobs.hsrl import SDS_HALO
                    obs_obj = HSRL(h5_file, Nav_only=True, verbose=False, SDS=SDS_HALO)
                    
                else:
                    if campaign.upper() in ['DISCOVERAQ', 'DISCOVER-AQ'] and plane.upper() == 'UC12':
                        # Custom mapping for DISCOVER-AQ HSRL-1 on the UC12 (omitted gps_date to bypass python3 bug)
                        SDS_DISCOVERAQ = {
                            'Nav_Data': ('gps_alt', 'gps_lat', 'gps_lon', 'gps_time'),
                            'Data_Products': ('Altitude', '532_ext', '532_bsc', '532_dep', 'Sa_532', '1064_bsc', '1064_ext')
                        }
                        obs_obj = HSRL(h5_file, Nav_only=True, verbose=False, SDS=SDS_DISCOVERAQ)
                        
                    elif campaign.upper() in ['DISCOVERAQ', 'DISCOVER-AQ'] and plane.upper() == 'B200':
                        # Custom mapping for DISCOVER-AQ HSRL-2 on the B200 (omitted Date to bypass python3 bug)
                        SDS_B200 = {
                            'ApplanixIMU': ('gps_alt', 'gps_lat', 'gps_lon', 'gps_time'),
                            'DataProducts': ('Altitude', '532_ext', '532_bsc', '532_dep', '532_Sa', '1064_bsc', '1064_ext')
                        }
                        obs_obj = HSRL(h5_file, Nav_only=True, verbose=False, SDS=SDS_B200)
                        
                    else:
                        from pyobs.hsrl import SDS_HSRL2
                        import copy
                        SDS_CUSTOM = copy.deepcopy(SDS_HSRL2)
                        if 'ApplanixIMU' in SDS_CUSTOM:
                            # This was a workaround for CAMP2EX
                            SDS_CUSTOM['Nav_Data'] = SDS_CUSTOM.pop('ApplanixIMU')
                        obs_obj = HSRL(h5_file, Nav_only=True, verbose=False, SDS=SDS_CUSTOM)
                    
                lon = obs_obj.lon[:].ravel()
                lat = obs_obj.lat[:].ravel()
                tyme = obs_obj.tyme.ravel()

            except Exception as e:
                print(f"Error reading {h5_file}: {e}")
                traceback.print_exc()
                continue

            if len(tyme) == 0:
                print(f"No valid time arrays found in {basename}. Skipping.")
                continue

            # Sort this specific file's trajectory
            sort_idx = np.argsort(tyme)
            lon = lon[sort_idx]
            lat = lat[sort_idx]
            tyme = tyme[sort_idx]
            
            # Calculate time step for this file
            median_dt = pd.Series(tyme).diff().dt.total_seconds().median()
            idt = int(60 / median_dt + 0.5) if pd.notna(median_dt) and median_dt > 0 else 1
            if idt < 1: 
                idt = 1
            lon, lat, tyme = lon[::idt], lat[::idt], tyme[::idt]

            # Extract suffix from filename (e.g., '_R2_L1_SUB' or '_R0_F1_')
            if date_str in basename:
                suffix = basename.split(date_str)[1].replace('.h5', '')
            else:
                suffix = "_R0" # Fallback

            # --- Sample MERRA-2 for this file ---
            if m2_files:
                out_m2 = f"{config['output_dir']}/{campaign}-MERRA2-aer-Nv-{plane}_Model_{date_str}{suffix}.nc"
                if not os.path.exists(out_m2):
                    print(f"  -> Sampling MERRA-2 (Using {len(m2_files)} daily files to bracket time)")
                    traj_m2 = TRAJECTORY(tyme, lon, lat, m2_files, verbose=True, chunks=chunks)
                    traj_ds_m2 = traj_m2.sample().compute()
                    traj_ds_m2.to_netcdf(out_m2, engine='netcdf4')
                else:
                    print(f"  -> MERRA-2 sampled file already exists: {out_m2}")
            else:
                print(f"  -> MERRA-2 files missing for bracketed dates.")
            
            # --- Sample MERRA-21C for this file ---
            if m21c_files:
                out_m21c = f"{config['output_dir']}/{campaign}-MERRA21C-aer-Nv-{plane}_Model_{date_str}{suffix}.nc"
                if not os.path.exists(out_m21c):
                    print(f"  -> Sampling MERRA-21C (Using {len(m21c_files)} 3-hourly files to bracket time)")
                    traj_m21c = TRAJECTORY(tyme, lon, lat, m21c_files, verbose=True, chunks=chunks)
                    traj_ds_m21c = traj_m21c.sample().compute()
                    traj_ds_m21c.to_netcdf(out_m21c, engine='netcdf4')
                else:
                    print(f"  -> MERRA-21C sampled file already exists: {out_m21c}")
            else:
                print(f"  -> Regridded MERRA-21C files missing for bracketed dates.")
