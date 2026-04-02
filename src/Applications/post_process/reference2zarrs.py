#!/usr/bin/env python3
"""
Convert an existing parquet reference store to per-variable Zarr stores.
Each variable is written independently — easy to parallelize and restart.
"""

import fsspec
import os, sys
import yaml
import argparse
import logging
import xarray as xr
from dask.distributed import Client, LocalCluster
import psutil
import time
from numcodecs import Blosc

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("config", help='configuration yaml file')
    parser.add_argument("--chunk_time", type=int, default=10,
                        help='number of time steps per chunk (default=10)')
    parser.add_argument("--n_workers", type=int, default=50,
                        help='number of dask workers (default=50)')
    parser.add_argument("--overwrite", action='store_true',
                        help='overwrite existing zarr files (default: skip)')
    parser.add_argument("--log_file", type=str, default=None,
                        help='optional log file path')
    parser.add_argument("--varname", type=str, default=None,
                        help='optional variable name')

    args = parser.parse_args()

    # setup logging
    fmt = '%(asctime)s - %(levelname)s - %(message)s'
    datefmt = '%Y-%m-%d %H:%M:%S'
    handlers = [logging.StreamHandler()]
    if args.log_file:
        handlers.append(logging.FileHandler(args.log_file))
    logging.basicConfig(level=logging.INFO, format=fmt,
                        datefmt=datefmt, handlers=handlers)
    logger = logging.getLogger(__name__)

    config = yaml.safe_load(open(args.config))
    com_outdir = config['combined_references_dir']

    ctl = config.get('reference_ctl')
    if not ctl:
        raise ValueError('No reference control files given. Cannot proceed.')

    # log memory info
    mem = psutil.virtual_memory()
    available_gb = mem.available / 1e9
    mem_per_worker = int(available_gb * 0.8 / args.n_workers)
    logger.info(f"Total node RAM:          {mem.total / 1e9:.1f} GB")
    logger.info(f"Available RAM:           {mem.available / 1e9:.1f} GB")
    logger.info(f"N workers:               {args.n_workers}")
    logger.info(f"Memory limit per worker: {mem_per_worker} GB")

    ref_parq = com_outdir + f'/reference_store_{ctl}.parq'

    # output directory for per-variable zarr stores
    zarr_dir = com_outdir + f'/zarr_store_{ctl}'
    os.makedirs(zarr_dir, exist_ok=True)

    if not os.path.exists(ref_parq):
        logger.error(f"Parquet reference store not found: "
                       f"{ref_parq} -- exiting")
        sys.exit(1)

    # =====================================================================
    # Step 1: Open the parquet reference store
    # =====================================================================
    logger.info(f"Opening parquet reference store: {ref_parq}")
    fs = fsspec.filesystem("reference", fo=ref_parq,
                           remote_protocol='file', lazy=True)
    ds = xr.open_dataset(fs.get_mapper(""), engine="zarr",
                         chunks={'time': args.chunk_time},
                         consolidated=False)

    logger.info(f"Dataset size:  {ds.nbytes / 1e9:.2f} GB")
    logger.info(f"Dataset dims:  {dict(ds.sizes)}")
    logger.info(f"Time range:    {ds.time.values[0]} to "
                f"{ds.time.values[-1]}")

    if args.varname is not None:
        Variables = [args.varname]
    else:
        Variables = list(ds.data_vars)
    n_vars = len(Variables)
    logger.info(f"Variables ({len(Variables)}): {Variables}")

    # =====================================================================
    # Step 2: Write one zarr per variable
    # =====================================================================
    local_dir = os.environ.get('LOCAL_TMPDIR', '/tmp')
    dask_kwargs = dict(n_workers=args.n_workers, threads_per_worker=1,
                       memory_limit=f'{mem_per_worker}GB',
                       local_directory=local_dir)

    with LocalCluster(**dask_kwargs) as cluster, Client(cluster) as client:

        for i, vn in enumerate(Variables):
            zarr_out = f'{zarr_dir}/{vn}.zarr'

            # skip if already done and not overwriting
            if os.path.exists(zarr_out) and not args.overwrite:
                logger.info(f"[{i+1}/{n_vars}] Skipping {vn} "
                            f"-- already exists: {zarr_out}")
                continue

            logger.info(f"[{i+1}/{n_vars}] Converting {vn}")
            logger.info(f"  Shape: {ds[vn].shape}")
            logger.info(f"  Size:  {ds[vn].nbytes / 1e9:.2f} GB")

            # set chunk sizes
            if 'lev' in ds[vn].dims:
                target_chunks = {'time': args.chunk_time, 'lev': -1,
                                 'lat': -1, 'lon': -1}
            else:
                target_chunks = {'time': args.chunk_time, 'lat': -1, 'lon': -1}


            da = ds[vn].chunk(target_chunks)
            # remove conflicting encoding attributes from coordinates
            for coord in da.coords:
                da[coord].attrs.pop('units', None)
                da[coord].attrs.pop('calendar', None)


            t0 = time.time()
            da.to_dataset(name=vn).to_zarr(
                zarr_out,
                mode='w',
                encoding={vn: {'compressor': Blosc(cname='lz4', clevel=5)}}
            )
            elapsed = time.time() - t0
            logger.info(f"  Written in {elapsed:.1f}s "
                        f"({ds[vn].nbytes / 1e9 / elapsed:.2f} GB/s)"
                        f" --> {zarr_out}")

    # =====================================================================
    # Summary
    # =====================================================================
    logger.info(f"Summary for {ctl}:")
    logger.info(f"  Per-variable zarr stores written to: {zarr_dir}/")
    for vn in Variables:
        zarr_out = f'{zarr_dir}/{vn}.zarr'
        if os.path.exists(zarr_out):
            vsize = sum(
                os.path.getsize(os.path.join(dp, f))
                for dp, _, files in os.walk(zarr_out)
                for f in files
            ) / 1e9
            logger.info(f"    {vn}.zarr  {vsize:.2f} GB")
        else:
            logger.warning(f"    {vn}.zarr  NOT FOUND -- conversion may "
                           f"have failed")

    logger.info(f"To use in sampling script:")
    logger.info(f"  import xarray as xr")
    logger.info(f"  ds = xr.open_mfdataset(")
    logger.info(f"      [f'{zarr_dir}/{{vn}}.zarr' for vn in Variables],")
    logger.info(f"      engine='zarr', chunks={{'time': {args.chunk_time}}}")
    logger.info(f"  )")
