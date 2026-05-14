#!/usr/bin/env python3
"""
    Samples MERRA-2 and MERRA-21C at 60 second intervals according to flight tracks in HSRL2 files.
"""
import sys
import os
import yaml
import numpy as np
from glob import glob
from datetime import datetime, timedelta
import pandas as pd

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
    obs_revision = config['obs_revision']

    kwargs = dict(n_workers=6, threads_per_worker=1, memory_limit='4GB')
    with LocalCluster(**kwargs) as cluster, Client(cluster) as client:
        
        h5_file = f"{config['obs_dir']}/{campaign.lower()}-{instrument}_{plane}_{date_str}_{obs_revision}.h5"
        
        if not os.path.exists(h5_file):
            print(f"Observation file missing: {h5_file}. Skipping sampling...")
            sys.exit(0)
            
        print(f"+++++++ Sampling on {h5_file}")
        
        try:
            if instrument.upper() == 'DIAL':
                obs_obj = DIAL(h5_file, verbose=False)
            else:
                obs_obj = HSRL(h5_file, Nav_only=True, verbose=False)
                
            lon = obs_obj.lon[:].ravel()
            lat = obs_obj.lat[:].ravel()
            tyme = obs_obj.tyme.ravel()
            
            median_dt = pd.Series(tyme).diff().dt.total_seconds().median()
            
            idt = int(60 / median_dt + 0.5)
            if idt < 1: 
                idt = 1
            lon, lat, tyme = lon[::idt], lat[::idt], tyme[::idt]

        except Exception as e:
            print(f"Error reading {h5_file}: {e}")
            sys.exit(0)

        chunks = {'time': 1, 'lev': -1, 'lat': -1, 'lon': -1}
        
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

        # --- Sample MERRA-2 ---
        if m2_files:
            out_m2 = f"{config['output_dir']}/{campaign}-MERRA2-aer-Nv-{plane}_Model_{date_str}_R0.nc"
            if not os.path.exists(out_m2):
                print(f"  -> Sampling MERRA-2 for {date_str} (Using {len(m2_files)} daily files to bracket time)")
                traj_m2 = TRAJECTORY(tyme, lon, lat, m2_files, verbose=True, chunks=chunks)
                traj_ds_m2 = traj_m2.sample().compute()
                traj_ds_m2.to_netcdf(out_m2, engine='netcdf4')
            else:
                print(f"  -> MERRA-2 sampled file already exists: {out_m2}")
        else:
            print(f"  -> MERRA-2 files missing for bracketed dates.")
        
        # --- Sample MERRA-21C ---
        if m21c_files:
            out_m21c = f"{config['output_dir']}/{campaign}-MERRA21C-aer-Nv-{plane}_Model_{date_str}_R0.nc"
            if not os.path.exists(out_m21c):
                print(f"  -> Sampling MERRA-21C for {date_str} (Using {len(m21c_files)} 3-hourly files to bracket time)")
                traj_m21c = TRAJECTORY(tyme, lon, lat, m21c_files, verbose=True, chunks=chunks)
                traj_ds_m21c = traj_m21c.sample().compute()
                traj_ds_m21c.to_netcdf(out_m21c, engine='netcdf4')
            else:
                print(f"  -> MERRA-21C sampled file already exists: {out_m21c}")
        else:
            print(f"  -> Regridded MERRA-21C files missing for bracketed dates.")
