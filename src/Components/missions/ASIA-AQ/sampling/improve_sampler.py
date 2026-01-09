#!/usr/bin/env python3
"""
    Sampler for IMPROVE surface PM2.5 obs for a modeling experiment
    This uses a parquet reference store to read the dataset
    Run ctl2kerchunk_run.j followed by kerchunk2parquet_run.j to create it
"""
import os, sys
import argparse
import yaml
from glob import glob
from datetime import datetime, timedelta
from pyobs.sampler import STATION
from pyobs.improve import SITE_MAP
from pyobs.xrctl import parse_ctl
from dask.distributed import performance_report
from dask.distributed import Client, LocalCluster
if __name__ == '__main__':


    parser = argparse.ArgumentParser()
    parser.add_argument("config",help='configuration yaml file')
    parser.add_argument("--profile",action="store_true")
    parser.add_argument("--n_workers",type=int,default=120,
                        help='number of pool workers to use (default=120)')
    parser.add_argument("--memory_limist",default='4BG',
                        help='memory limit per worker (default=4GB)')

    args = parser.parse_args()

    config = yaml.safe_load(open(args.config))

    # get IMPROVE site locations
    site_path = config['improve_site_map']
    sites = SITE_MAP(site_path)

    # create output directory
    outdir = config['sampled_outdir'] + '/improve'
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    # get reference parquet location
    refdir = config['combined_references_dir']


    kwargs = dict(n_workers=args.n_workers, threads_per_worker=1, memory_limit=args.memory_limit)
    with LocalCluster(**kwargs) as cluster, Client(cluster) as client:
        print(f"+++++++ IMPROVE Sampling on {site_path}")
        station,lat,lon = sites.df['SiteCode'].values, sites.df['Longitude'].values, sites.df['Latitude'].values

        # Sample Aerosol Collection
        # --------------------------------------
        ctls = [config['model_aer_ctl']] + config['model_other_ctl']
        ctls = [x for x in ctls if x is not None]
#        chunks = {'time':1, 'lev':-1, 'lat':-1, 'lon': -1}
        chunks = 'auto'
        for ctl in ctls:
            # get reference parquet
            outparq = refdir + f'/reference_store_{ctl}.parq'
                 
            print(f" Sampling on reference {outparq}")                
            stn = STATION(station,lon,lat,outparq,verbose=True,chunks=chunks,engine='h5netcdf')
                
            stn_ds = stn.sample()

            # write out the native sampled model fields
            outFile = f'{outdir}/improve.{ctl}.nc4'
            if args.profile:
                with performance_report(filename="dask-report.html"):
                    stn_ds = stn_ds.compute()
            else:
                stn_ds = stn_ds.compute()
            stn_ds.to_netcdf(outFile,format='NETCDF4_CLASSIC')
            print(f"Successfully wrote {outFile} for {ctl} control file")

            client.wait_for_workers(args.n_workers)
