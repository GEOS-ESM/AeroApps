#!/usr/bin/env python3
"""
Convert an existing parquet reference store to a consolidated Zarr store.
"""

import fsspec
import os, sys
import yaml
import argparse
import xarray as xr
from dask.distributed import Client, LocalCluster
import psutil

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("config", help='configuration yaml file')
    parser.add_argument("--chunk_time", type=int, default=10,
                        help='number of time steps per chunk (default=10)')
    parser.add_argument("--lat_min", type=float, default=None,
                        help='minimum latitude for spatial subset (optional)')
    parser.add_argument("--lat_max", type=float, default=None,
                        help='maximum latitude for spatial subset (optional)')
    parser.add_argument("--lon_min", type=float, default=None,
                        help='minimum longitude for spatial subset (optional)')
    parser.add_argument("--lon_max", type=float, default=None,
                        help='maximum longitude for spatial subset (optional)')
    parser.add_argument("--n_workers", type=int, default=50,
                        help='number of dask workers (default=50)')

    args = parser.parse_args()

    config = yaml.safe_load(open(args.config))

    com_outdir = config['combined_references_dir']

    ctls = config['reference_ctl_list']
    if len(ctls) == 0:
        raise ValueError('No reference control files given. Cannot proceed.')

    # log memory info
    mem = psutil.virtual_memory()
    print(f"Total node RAM:  {mem.total / 1e9:.1f} GB")
    print(f"Available RAM:   {mem.available / 1e9:.1f} GB")
    print(f"N workers:       {args.n_workers}")

    # set memory limit
    available_gb = mem.available / 1e9
    mem_per_worker = int(available_gb * 0.8 / args.n_workers)
    print(f"Memory limit per worker: {mem_per_worker} GB")

    for ctl in ctls:
        ref_parq = com_outdir + f'/reference_store_{ctl}.parq'
        zarr_out = com_outdir + f'/zarr_store_{ctl}.zarr'

        if not os.path.exists(ref_parq):
            print(f"WARNING: Parquet reference store not found: {ref_parq} -- skipping")
            continue

        print(f"")
        print(f"Converting {ref_parq}")
        print(f"       --> {zarr_out}")

        # =====================================================================
        # Step 1: Open the parquet reference store
        # =====================================================================
        print(f"Opening parquet reference store...")
        fs = fsspec.filesystem("reference", fo=ref_parq,
                               remote_protocol='file', lazy=True)
        ds = xr.open_dataset(fs.get_mapper(""), engine="zarr",
                             chunks={'time': args.chunk_time},
                             consolidated=False)

        print(f"Dataset size:  {ds.nbytes / 1e9:.2f} GB")
        print(f"Dataset dims:  {dict(ds.sizes)}")
        print(f"Variables:     {len(ds.data_vars)}")
        print(f"Time range:    {ds.time.values[0]} to {ds.time.values[-1]}")

        # =====================================================================
        # Step 2: Optional spatial subset
        # =====================================================================
        if all(v is not None for v in [args.lat_min, args.lat_max,
                                        args.lon_min, args.lon_max]):
            print(f"Subsetting to lat [{args.lat_min}, {args.lat_max}], "
                  f"lon [{args.lon_min}, {args.lon_max}]")
            ds = ds.sel(lat=slice(args.lat_min, args.lat_max),
                        lon=slice(args.lon_min, args.lon_max))
            print(f"Subsetted size: {ds.nbytes / 1e9:.2f} GB")
            print(f"Subsetted dims: {dict(ds.sizes)}")

        # =====================================================================
        # Step 3: Rechunk to good chunk sizes for zarr
        # =====================================================================
        if 'lev' in ds.dims:
            target_chunks = {'time': args.chunk_time, 'lev': -1, 'lat': -1, 'lon': -1}
        else:
            target_chunks = {'time': args.chunk_time, 'lat': -1, 'lon': -1}

        ds = ds.chunk(target_chunks)
        chunk_summary = {dim: c[0] for dim, c in ds.chunks.items()}
        print(f"Target chunks: {chunk_summary}")

        # estimate single chunk size
        vn0 = list(ds.data_vars)[0]
        chunk_size_mb = ds[vn0].isel(time=0).nbytes / 1e6
        print(f"Single time step size: {chunk_size_mb:.2f} MB")

        # =====================================================================
        # Step 4: Write to zarr using Dask
        # =====================================================================
        local_dir = os.environ.get('LOCAL_TMPDIR', '/tmp')
        dask_kwargs = dict(n_workers=args.n_workers, threads_per_worker=1,
                           memory_limit=f'{mem_per_worker}GB',
                           local_directory=local_dir)

        print(f"Writing zarr store...")
        with LocalCluster(**dask_kwargs) as cluster, Client(cluster) as client:
            # no compression for maximum read speed
            encoding = {vn: {'compressor': None} for vn in ds.data_vars}
            ds.to_zarr(zarr_out, mode='w', encoding=encoding)
            print(f"Successfully wrote zarr store: {zarr_out}")

        # =====================================================================
        # Summary
        # =====================================================================
        zarr_size = sum(
            os.path.getsize(os.path.join(dirpath, f))
            for dirpath, _, filenames in os.walk(zarr_out)
            for f in filenames
        ) / 1e9

        print(f"")
        print(f"Summary:")
        print(f"  Input parquet:  {ref_parq}")
        print(f"  Output zarr:    {zarr_out}")
        print(f"  Zarr size:      {zarr_size:.2f} GB")
        print(f"")
        print(f"To use in sampling script:")
        print(f"  ds = xr.open_zarr('{zarr_out}', chunks={{'time': {args.chunk_time}}})")
