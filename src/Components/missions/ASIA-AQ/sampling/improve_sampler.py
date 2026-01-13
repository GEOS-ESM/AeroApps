#!/usr/bin/env python3
"""
    Sampler for IMPROVE surface PM2.5 obs for a modeling experiment
    This uses a parquet reference store to read the dataset
    Run ctl2kerchunk_run.j followed by kerchunk2parquet_run.j to create it
"""
import os, sys
import argparse
import yaml
from pyobs.improve import SITE_MAP
from dask.distributed import Client, LocalCluster
import fsspec
import xarray as xr
import xesmf as xe
import dask

# Disable worker heartbeats to prevent shutdown race conditions
dask.config.set({"distributed.scheduler.worker-ttl": None})

if __name__ == '__main__':


    parser = argparse.ArgumentParser()
    parser.add_argument("config",help='configuration yaml file')
    parser.add_argument("--n_workers",type=int,default=50,
                        help='number of pool workers to use (default=50)')
    parser.add_argument("--memory_limit",default='9GB',
                        help='memory limit per worker (default=9GB)')

    args = parser.parse_args()

    config = yaml.safe_load(open(args.config))

    # get IMPROVE site locations
    site_path = config['improve_site_map']
    sites = SITE_MAP(site_path)
    station,lon,lat = sites.df['SiteCode'].values, sites.df['Longitude'].values, sites.df['Latitude'].values
    lon = xr.DataArray(lon, dims='station')
    lat = xr.DataArray(lat, dims='station')
    ds_loc = xr.Dataset({"lon": lon, "lat": lat})

    # create output directory
    outdir = config['sampled_outdir'] + '/improve'
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    # get reference parquet location
    refdir = config['combined_references_dir']


    kwargs = dict(n_workers=args.n_workers, threads_per_worker=1, memory_limit=args.memory_limit)
    with LocalCluster(**kwargs) as cluster, Client(cluster) as client:
        print(f"+++++++ IMPROVE Sampling on {site_path}")

        # Sample Aerosol Collection
        # --------------------------------------
        ctls = [config['model_aer_ctl']] + config['model_other_ctl']
        ctls = [x for x in ctls if x is not None]

#        chunks = {'time':10, 'lev':-1, 'lat':-1, 'lon': -1}
#        chunks = {'time':10}
        chunks = 'auto'

        for ctl in ctls:
            # get reference parquet
            path_to_parq = refdir + f'/reference_store_{ctl}.parq'
            print(f" Sampling on reference {path_to_parq}")
            fs = fsspec.filesystem("reference", fo=path_to_parq, remote_protocol='file', lazy=True)
            ds = xr.open_dataset(fs.get_mapper(""), engine="zarr", chunks=chunks, consolidated=False)
            Variables = list(ds.data_vars)
            regridder = xe.Regridder(ds, ds_loc, "bilinear", locstream_out=True)
            
            outFile = f'{outdir}/improve.{ctl}.nc4'
            first = True
            for vn in Variables:
                print(f'Sampling {vn}')
                sampled = {}
                stn_ds = regridder(ds[vn])
                sampled[vn] = stn_ds.compute()
                stn_ds = xr.Dataset(sampled).assign_coords({'station': station})
                # write out the native sampled model fields
                if first:
                    stn_ds.to_netcdf(outFile,format='NETCDF4_CLASSIC')
                    first = False
                else:
                    stn_ds.to_netcdf(outFile,format='NETCDF4_CLASSIC',mode="a")    
                print(f"Successfully wrote {vn} to {outFile}")


            print(f"Successfully wrote {outFile} for {ctl} control file")

            client.wait_for_workers(args.n_workers)
