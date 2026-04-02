#!/usr/bin/env python3
"""
    Sampler for AMON surface NH4 and NH3 obs for a modeling experiment
    This uses a parquet reference store to read the dataset
    Run ctl2kerchunk_run.j followed by kerchunk2parquet_run.j to create it
"""
import os, sys
import argparse
import yaml
from dask.distributed import Client, LocalCluster
import fsspec
import xarray as xr
import xesmf as xe
import dask
import pandas as pd
import logging
from distributed import wait
import gc
import time
import psutil
import json

# Disable worker heartbeats to prevent shutdown race conditions
dask.config.set({"distributed.scheduler.worker-ttl": None})

if __name__ == '__main__':


    parser = argparse.ArgumentParser()
    parser.add_argument("config",help='configuration yaml file')
    parser.add_argument("--n_workers",type=int,default=50,
                        help='number of pool workers to use (default=50)')

    args = parser.parse_args()

    config = yaml.safe_load(open(args.config))

    # setup logger
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    logger = logging.getLogger(__name__)

    # get AMoN site locations
    site_path = config['amon_site_map']
    sites = pd.read_table(site_path, sep=',')
    sites = sites.dropna(subset=['latitude', 'longitude'])  # drop sites where lat or lon are nan
    station,lon,lat = sites['siteId'].values, sites['longitude'].values, sites['latitude'].values
    lon = xr.DataArray(lon, dims='station')
    lat = xr.DataArray(lat, dims='station')
    ds_loc = xr.Dataset({"lon": lon, "lat": lat})

    # Compute station bounding box 
    lat_min = float(ds_loc.lat.values.min()) - 1
    lat_max = float(ds_loc.lat.values.max()) + 1
    lon_min = float(ds_loc.lon.values.min()) - 1
    lon_max = float(ds_loc.lon.values.max()) + 1

    # Log bounding box
    logger.info(f"Station bounding box: lat [{lat_min:.1f}, {lat_max:.1f}], "
                f"lon [{lon_min:.1f}, {lon_max:.1f}]")

    # create output directory
    outdir = config['sampled_outdir'] + '/amon'
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    # get reference parquet location
    refdir = config['combined_references_dir']

    # Check total available memory on node
    mem = psutil.virtual_memory()
    logger.info(f"Total node RAM: {mem.total / 1e9:.1f} GB")
    logger.info(f"Available RAM: {mem.available / 1e9:.1f} GB")
    logger.info(f"N workers: {args.n_workers}")
    logger.info(f"Memory per worker: {mem.total / args.n_workers / 1e9:.1f} GB")

    available_gb = mem.available / 1e9

    # set memory limits
    mem_per_worker = int(available_gb * 0.8 / args.n_workers)
    local_dir = os.environ.get('LOCAL_TMPDIR', '/tmp')
    kwargs = dict(n_workers=args.n_workers, threads_per_worker=1, 
                  memory_limit=f'{mem_per_worker}GB',local_directory=local_dir)

    chunks = {'time':1, 'lev':-1, 'lat':-1, 'lon': -1}

    logger.info(f"Config: n_workers={args.n_workers}, "
            f"mem_limit={mem_per_worker}GB, "
            f"chunk_time={chunks['time']}")

    with LocalCluster(**kwargs) as cluster, Client(cluster) as client:
        logger.info(f"+++++++ AMoN Sampling on {site_path}")

        # Sample Aerosol Collection
        # --------------------------------------
        ctls = [config['model_aer_ctl']] + config['model_other_ctl']
        ctls = [x for x in ctls if x is not None]

        for ctl in ctls:
            # get reference parquet
            path_to_parq = refdir + f'/reference_store_{ctl}.parq'

            # Read the parquet reference store to see the source files
            ref_df = pd.read_parquet(path_to_parq)
            logger.info(f"Reference store shape: {ref_df.shape}")
            logger.info(f"Reference store columns: {ref_df.columns.tolist()}")
            logger.info(f"\n{ref_df.head(10)}")

            # How many unique files are being referenced?
            if 'path' in ref_df.columns:
                unique_files = ref_df['path'].nunique()
                logger.info(f"Number of unique source files: {unique_files}")

            logger.info(f" Sampling on reference {path_to_parq}")
            fs = fsspec.filesystem("reference", fo=path_to_parq, remote_protocol='file', lazy=True)
            ds = xr.open_dataset(fs.get_mapper(""), engine="zarr", chunks=chunks, consolidated=False)
            ds = ds.unify_chunks()

            # SPATIAL SUBSET — reduces data reads dramatically
            logger.info(f"Full dataset size: {ds.nbytes / 1e9:.2f} GB")
            ds_sub = ds.sel(
                lat=slice(lat_min, lat_max),
                lon=slice(lon_min, lon_max)
            )
            logger.info(f"Subsetted dataset size: {ds_sub.nbytes / 1e9:.2f} GB "
                        f"({ds_sub.nbytes / ds.nbytes * 100:.1f}% of original)")
            logger.info(f"Subsetted dims: {dict(ds_sub.sizes)}")

            Variables = list(ds_sub.data_vars)

            logger.info(f"Dataset size: {ds_sub.nbytes / 1e9:.2f} GB")
            logger.info(f"Dataset dims: {dict(ds_sub.sizes)}")
            chunk_summary = {dim: chunks[0] for dim, chunks in ds_sub.chunks.items()}
            logger.info(f"Chunks: {chunk_summary}")
            logger.info(f"Number of variables: {len(ds_sub.data_vars)}")


            vn = Variables[0]
            logger.info(f"Single variable size: {ds_sub[vn].nbytes / 1e9:.2f} GB")
            logger.info(f"Single variable shape: {ds_sub[vn].shape}")
            logger.info(f"Single time step size: {ds_sub[vn].isel(time=0).nbytes / 1e6:.2f} MB")            

            t0 = time.time()
            regridder = xe.Regridder(ds_sub, ds_loc, "bilinear", locstream_out=True)
            logging.info(f'Computed regidder in {time.time()-t0:.2f}s')

            outFile = f'{outdir}/amon.{ctl}.nc4'
            First = True
            # Process variables one at a time
            for vn in Variables:
                logger.info(f'Sampling {vn} - shape {ds_sub[vn].shape}')

                t0 = time.time()
                # Persist variable in distributed memory before computing
                # This reads the data once and keeps it in worker RAM
                da_persisted = client.persist(ds_sub[vn])

                # Wait for all workers to finish loading
                wait(da_persisted)
                logger.info(f"Persisted {vn} in {time.time()-t0:.2f}s")


                t0 = time.time()
                stn_da = regridder(da_persisted)
                logger.info(f"  Output shape after regrid: {stn_da.shape}")
                logger.info(f"  Output size: {stn_da.nbytes / 1e6:.2f} MB")  # should be tiny
 
                stn_da = stn_da.compute()
                logger.info(f"Computed {vn} in {time.time()-t0:.2f}s")

                t0 = time.time()
                stn_ds = xr.Dataset({vn: stn_da}).assign_coords({'station': station})
                if First:
                    stn_ds.to_netcdf(outFile, format='NETCDF4_CLASSIC')
                    first = False
                else:
                    stn_ds.to_netcdf(outFile, format='NETCDF4_CLASSIC', mode='a')
                logger.info(f"Wrote {vn} in {time.time()-t0:.2f}s")

            logger.info(f"Successfully wrote {outFile} for {ctl} control file")

            # give workers time to shut down cleanly
            client.run(gc.collect)
            client.shutdown()   # explicitly shut down workers before cluster exits
