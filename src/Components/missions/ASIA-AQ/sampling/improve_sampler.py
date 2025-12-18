#!/usr/bin/env python3
"""
    Sampler for IMPROVE surface PM2.5 obs for a modeling experiment
"""
import os, sys
import argparse
import yaml
from glob import glob
from datetime import datetime
from pyobs.sampler import STATION
from pyobs.improve import SITE_MAP
from pyobs.xrctl import parse_ctl
from dask.distributed import performance_report
from dask.distributed import Client, LocalCluster
if __name__ == '__main__':


    parser = argparse.ArgumentParser()
    parser.add_argument("config",help='configuration yaml file')
    parser.add_argument("--profile",action="store_true")

    args = parser.parse_args()

    config = yaml.safe_load(open(args.config))

    # create output directory
    if not os.path.exists(config['sampled_outdir']):
        os.makedirs(config['sampled_outdir'])

    # get IMPROVE site locations 
    site_path = config['improve_site_map']
    sites = SITE_MAP(site_path)

    # get years to sample
    years = config['improve_years']

    kwargs = dict(n_workers=120, threads_per_worker=1, memory_limit='4GB')
    with LocalCluster(**kwargs) as cluster, Client(cluster) as cluster:
        for year in years:
            print(f"+++++++ IMPROVE Sampling on {site_path}")
            time_range = [datetime(year,1,1),datetime(year,12,31)]
            station,lon,lat = sites.df['SiteCode'],sites.df['Longitude'],sites.df['Latitude']

            # Sample Aerosol Collection
            # --------------------------------------
            ctl = config['model_aer_ctl'] 
            chunks = {'time':1, 'lev':-1, 'lat':-1, 'lon': -1}
            # check that all files exist
            paths = parse_ctl(ctl,time_range)
            exists = True
            for p in paths:
                if not os.path.exists(p): exists = False
            
            if exists:
                stn = STATION(station,lon,lat,ctl,time_range=time_range,verbose=True,chunks=chunks,engine='h5netcdf')
                
                stn_ds = stn.sample()
                if args.profile:
                    with performance_report(filename="dask-report.html"):
                        stn_ds = stn_ds.compute()
                else:
                    stn_ds = stn_ds.compute()

                # write out the native sampled model fields
                outFile = config['sampled_outdir'] + '/improve.'  + ctl + '.nc4'
                stn_ds.to_netcdf(outFile,engine='netcdf4')

                # If other model collections are provided
                # sample those and write out
                # -----------------------------
                for ctl in config['model_other_ctl']:
                    if ctl is not None:
                        try:
                            # check that all files exist
                            paths = parse_ctl(ctl,time_range)
                            exists = True
                            for p in paths:
                                if not os.path.exists(p): exists = False

                            if exists:
                                stn = STATION(station,lon,lat,ctl,time_range=time_range,verbose=True,chunks=chunks,engine='h5netcdf')
                            else:
                                print(f"Model input files missing for {year}")
                                print("Skipping...")
                                continue
                        except:
                            print(f"Error reading model for {year}")
                            print("Skipping....")
                            continue

                        stn_ds = stn.sample()
                        stn_ds = stn_ds.compute()

                        # write out the native sampled model fields
                        outFile = config['sampled_outdir'] + '/improve.' + ctl + '.nc4'
                        stn_ds.to_netcdf(outFile,engine='netcdf4')

            else:
                print(f"Model input files missing for {year}")
                print("Skipping...")
                continue
