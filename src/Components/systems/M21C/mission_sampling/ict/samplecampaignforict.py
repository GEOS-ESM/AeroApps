#!/usr/bin/env python3
"""
    Generic GEOS Sampler for Airborne Field Campaigns 
"""
import os, sys
import argparse
import yaml
from glob import glob
import numpy as np
import pandas as pd
from pyobs.sampler import TRAJECTORY
from pyobs.icartt import ICARTT
from pyobs.xrctl import parse_ctl
from dask.distributed import performance_report
from dask.distributed import Client, LocalCluster

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("config", help='configuration yaml file')
    parser.add_argument("--profile", action="store_true")

    args = parser.parse_args()

    config = yaml.safe_load(open(args.config))

    # create output directory
    if not os.path.exists(config['sampled_outdir']):
        os.makedirs(config['sampled_outdir'])

    # get flight dates from config
    flight_dates = config.get('flight_dates', [])

    # get all ict files
    all_ictFiles = sorted(glob(config['mergefiles']+'/*ict'))

    # filter ict files based on flight_dates in the config
    if flight_dates:
        ictFiles = [f for f in all_ictFiles if any(date in os.path.basename(f) for date in flight_dates)]
        print(f"Found {len(ictFiles)} files matching the configured flight dates.")
    else:
        ictFiles = all_ictFiles
        print(f"No flight dates specified in config. Processing all {len(ictFiles)} files.")

    kwargs = dict(n_workers=40, threads_per_worker=1, memory_limit='4GB')
    
    with LocalCluster(**kwargs) as cluster, Client(cluster) as client:
        for ict in ictFiles:
            print(f"Working on {ict}")
            m = ICARTT(ict)
            lon, lat, tyme = getattr(m, config['lon']), getattr(m, config['lat']), m.Nav['Time']
            lon = np.where(lon > 180, lon - 360, lon)


            # Sample Aerosol Collection
            # --------------------------------------
            ctl = config['model_aer_ctl'] 
            chunks = {'time':1, 'lev':-1, 'lat':-1, 'lon': -1}
            
            time_range = [tyme.min(), tyme.max()]
            
            paths = parse_ctl(ctl, time_range)
            exists = len(paths) > 0 
            
            for p in paths:
                if not os.path.exists(p): 
                    print(f"  -> WARNING: Path does NOT exist: {p}")
                    exists = False
            
            if exists:
                traj = TRAJECTORY(
                    tyme, lon, lat, ctl, 
                    verbose=True, 
                    chunks=chunks, 
                    engine='h5netcdf', 
                    combine='nested', 
                    concat_dim='time',
                    data_vars='minimal',
                    join='override'
                )
                
                traj_ds = traj.sample()
                if args.profile:
                    with performance_report(filename="dask-report.html"):
                        traj_ds = traj_ds.compute()
                else:
                    traj_ds = traj_ds.compute()

                var_names = list(traj_ds.data_vars.keys())
                if var_names:
                    sample_var = var_names[0]

                # write out the native sampled model fields
                outFile = config['sampled_outdir'] + '/' + os.path.basename(ict)[:-3] + ctl + '.nc4'
                traj_ds.to_netcdf(outFile, engine='netcdf4')

                # If other model collections are provided
                # sample those and write out
                # -----------------------------
                if config['model_other_ctl'] is not None:
                    for ctl in config['model_other_ctl']:
                        if ctl is not None:
                            try:
                                time_range = [tyme.min(), tyme.max()]
                                paths = parse_ctl(ctl, time_range)
                                exists = len(paths) > 0
                                for p in paths:
                                    if not os.path.exists(p): exists = False

                                if exists:
                                    traj = TRAJECTORY(tyme, lon, lat, ctl, verbose=True, chunks=chunks, engine='h5netcdf', combine='nested', concat_dim='time')
                                else:
                                    print(f"Model input files missing for {ict} (ctl: {ctl})")
                                    print("Skipping...")
                                    continue
                            except Exception as e:
                                print(f"Error reading model for {ict} (ctl: {ctl}): {e}")
                                print("Skipping....")
                                continue

                            traj_ds = traj.sample()
                            traj_ds = traj_ds.compute()

                            outFile = config['sampled_outdir'] + '/' + os.path.basename(ict)[:-3] + ctl + '.nc4'
                            traj_ds.to_netcdf(outFile, engine='netcdf4')

            else:
                print(f"\nERROR: Model input files missing or parse_ctl returned empty list for {ict}")
                print("Skipping...")
                continue
