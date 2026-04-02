#!/usr/bin/env python3
"""
Create a consolidated Zarr store for a GEOS collection.
Maintains both:
  - parquet reference store (for backup/reconstruction)
  - consolidated zarr store (for fast sampling)

This script is still experimental.  It might kill memory as written.
"""

import fsspec
from fsspec.implementations.reference import LazyReferenceMapper
import os, sys
from pyobs.xrctl import parse_ctl
import yaml
import argparse
import numpy as np
from datetime import datetime
import kerchunk.hdf
import kerchunk.combine
from dask.distributed import Client, LocalCluster
import xarray as xr
import multiprocessing

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("config", help='configuration yaml file')
    parser.add_argument("--append", action='store_true',
                        help='append new timesteps to existing parquet and zarr stores')
    parser.add_argument("--n_workers", type=int, default=50,
                        help='number of pool workers to use (default=50)')
    parser.add_argument("--max_batches", type=int, default=100,
                        help='max number of batches (default=100)')
    parser.add_argument("--chunk_time", type=int, default=1,
                        help='number of time steps per chunk (default=1)')
    parser.add_argument("--lat_min", type=float, default=None,
                        help='minimum latitude for spatial subset (optional)')
    parser.add_argument("--lat_max", type=float, default=None,
                        help='maximum latitude for spatial subset (optional)')
    parser.add_argument("--lon_min", type=float, default=None,
                        help='minimum longitude for spatial subset (optional)')
    parser.add_argument("--lon_max", type=float, default=None,
                        help='maximum longitude for spatial subset (optional)')

    args = parser.parse_args()

    config = yaml.safe_load(open(args.config))

    ctls = config['reference_ctl_list']
    if len(ctls) == 0:
        raise ValueError('No reference control files given. Cannot proceed.')

    com_outdir = config['combined_references_dir']
    os.makedirs(com_outdir, exist_ok=True)

    fs_local = fsspec.filesystem("file")

    for ctl in ctls:
        # get the list of files on disk
        paths = parse_ctl(ctl)
        path_list = []
        for p in paths:
            if os.path.exists(p): path_list.append(p)
        path_list = np.array(path_list)

        if len(path_list) == 0:
            raise ValueError(f'No {ctl} files found on disk. Cannot proceed.')

        # both output stores
        ref_parq = com_outdir + f'/reference_store_{ctl}.parq'
        zarr_out = com_outdir + f'/zarr_store_{ctl}.zarr'

        # if appending, figure out which files are new
        if args.append:
            if not os.path.exists(ref_parq):
                raise ValueError(f'Parquet reference store {ref_parq} does not exist. '
                                 f'Run without --append first.')
            if not os.path.exists(zarr_out):
                raise ValueError(f'Zarr store {zarr_out} does not exist. '
                                 f'Run without --append first.')

            # use parquet to find last time — it is the source of truth
            fs = fsspec.filesystem("reference", fo=ref_parq,
                                   remote_protocol='file', lazy=True)
            ds_existing = xr.open_dataset(fs.get_mapper(), engine="zarr",
                                          chunks='auto',
                                          backend_kwargs={"consolidated": False})
            last_time = ds_existing.time[-1].values
            last_time = last_time.astype('M8[us]').astype('O')
            ds_existing.close()
            print(f"Last time in reference store: {last_time}")

            # keep only files after last time
            date_list = []
            for p in path_list:
                date_list.append(datetime.strptime(
                    os.path.basename(p).split('.')[-2], '%Y%m%d_%H%Mz'))
            date_list = np.array(date_list)
            keep = date_list > last_time
            path_list = path_list[keep]

            if len(path_list) == 0:
                print(f'No new {ctl} files found to append. Nothing to be done.')
                sys.exit(2)

            print(f"Found {len(path_list)} new files to append")

        # define function to convert HDF5 to kerchunk reference
        def process_single_file(path):
            try:
                with fsspec.open(path, mode='rb') as f:
                    scanner = kerchunk.hdf.SingleHdf5ToZarr(
                        f, path, inline_threshold=0).translate()
                return scanner
            except Exception as e:
                raise RuntimeError(f"Error processing {path}: {str(e)}")

        # Define the dimension you want to concatenate along (e.g., 'time')
        concat_dims = ['time']

        # List dimensions that are identical across files and should NOT be concatenated
        # figure out if there is a lev dimension
        with open(ctl) as f:
            content = f.read()
        if 'zdef' in content:
            identical_dims = ['lat', 'lon', 'lev']
        else:
            identical_dims = ['lat', 'lon']

        # 'cf:time' decodes internal 'time' variable using CF-conventions
        coo_map = {"time": "cf:time"}
        # 'M8[ns]' forces the result to NumPy datetime64 nanoseconds
        coo_dtypes = {"time": "M8[ns]"}

        # split files into batches
        path_split = np.array_split(
            path_list, max([1, len(path_list) // args.n_workers]))
        if len(path_split) > args.max_batches:
            print(f'{len(path_split)} batches found. Limiting to {args.max_batches}')
            path_split = path_split[:args.max_batches]
        path_split = [list(x) for x in path_split]

        # =====================================================================
        # Step 1: Build or append to the parquet reference store
        #         AND collect references for zarr write
        # =====================================================================
        # We accumulate all new single_refs so we can write to zarr
        # without re-reading from the parquet store
        all_new_refs = []

        if not args.append:
            # create the parquet reference store from the first batch
            out = LazyReferenceMapper.create(root=ref_parq, fs=fs_local,
                                             record_size=100000)
            first_batch = path_split.pop(0)
            n_workers = min([args.n_workers, len(first_batch)])

            with multiprocessing.Pool(n_workers) as pool:
                single_refs = pool.map(process_single_file, first_batch)

            mzz = kerchunk.combine.MultiZarrToZarr(
                path=single_refs,
                remote_protocol='file',
                concat_dims=concat_dims,
                identical_dims=identical_dims,
                coo_map=coo_map,
                coo_dtypes=coo_dtypes,
                out=out
            ).translate()
            out.flush()
            all_new_refs.extend(single_refs)
            del single_refs, mzz
            print(f"Successfully created parquet reference store: {ref_parq}")

        # process remaining batches — append to parquet and collect refs
        i = 0
        nbatches = len(path_split)
        with multiprocessing.Pool(args.n_workers) as pool:
            for current_batch in path_split:
                single_refs = pool.map(process_single_file, current_batch)
                print(f'Batch {i+1} of {nbatches} single references created')

                out = LazyReferenceMapper(root=ref_parq, fs=fs_local)
                mzz = kerchunk.combine.MultiZarrToZarr.append(
                    path=single_refs,
                    original_refs=out,
                    remote_protocol='file',
                    concat_dims=concat_dims,
                    identical_dims=identical_dims,
                    coo_map=coo_map,
                    coo_dtypes=coo_dtypes
                ).translate()
                out.flush()
                all_new_refs.extend(single_refs)

                print(f"Batch {i+1} of {nbatches} appended to parquet: {ref_parq}")
                i += 1
                del single_refs, mzz

        print(f"Successfully completed parquet reference store: {ref_parq}")

        # =====================================================================
        # Step 2: Open the new references directly and write/append to zarr
        #         Uses all_new_refs collected above — no need to re-read parquet
        # =====================================================================
        print(f"Opening new references and writing to zarr store: {zarr_out}")

        # combine all new refs into a single virtual dataset
        combined = kerchunk.combine.MultiZarrToZarr(
            path=all_new_refs,
            remote_protocol='file',
            concat_dims=concat_dims,
            identical_dims=identical_dims,
            coo_map=coo_map,
            coo_dtypes=coo_dtypes
        ).translate()

        fs_new = fsspec.filesystem("reference", fo=combined,
                                   remote_protocol='file')
        ds_new = xr.open_dataset(fs_new.get_mapper(""), engine="zarr",
                                 chunks={'time': args.chunk_time},
                                 consolidated=False)

        # spatial subset — only on first creation
        # on append the zarr already has the right spatial extent
        if not args.append and all(v is not None for v in [args.lat_min,
                                                            args.lat_max,
                                                            args.lon_min,
                                                            args.lon_max]):
            print(f"Subsetting to lat [{args.lat_min}, {args.lat_max}], "
                  f"lon [{args.lon_min}, {args.lon_max}]")
            ds_new = ds_new.sel(lat=slice(args.lat_min, args.lat_max),
                                lon=slice(args.lon_min, args.lon_max))
            print(f"Subsetted size: {ds_new.nbytes / 1e9:.2f} GB")

        # set chunk sizes
        if 'lev' in ds_new.dims:
            target_chunks = {'time': args.chunk_time, 'lev': -1, 'lat': -1, 'lon': -1}
        else:
            target_chunks = {'time': args.chunk_time, 'lat': -1, 'lon': -1}

        ds_new = ds_new.chunk(target_chunks)
        chunk_summary = {dim: c[0] for dim, c in ds_new.chunks.items()}
        print(f"Chunks: {chunk_summary}")
        print(f"Size to write: {ds_new.nbytes / 1e9:.2f} GB")

        # =====================================================================
        # Step 3: Write to zarr using Dask
        # =====================================================================
        local_dir = os.environ.get('LOCAL_TMPDIR', '/tmp')
        dask_kwargs = dict(n_workers=args.n_workers, threads_per_worker=1,
                           local_directory=local_dir)

        with LocalCluster(**dask_kwargs) as cluster, Client(cluster) as client:
            if args.append:
                print(f"Appending new timesteps to zarr store: {zarr_out}")
                ds_new.to_zarr(zarr_out, mode='a', append_dim='time')
            else:
                encoding = {vn: {'compressor': None} for vn in ds_new.data_vars}
                print(f"Writing new zarr store: {zarr_out}")
                ds_new.to_zarr(zarr_out, mode='w', encoding=encoding)

            print(f"Successfully wrote zarr store: {zarr_out}")

        # =====================================================================
        # Summary
        # =====================================================================
        print(f"")
        print(f"Summary:")
        print(f"  Parquet reference store (for backup/append): {ref_parq}")
        print(f"  Zarr store (for fast sampling):              {zarr_out}")
        print(f"")
        print(f"To append new files:")
        print(f"  python {sys.argv[0]} {args.config} --append")
        print(f"")
        print(f"To use in sampling script:")
        print(f"  ds = xr.open_zarr('{zarr_out}', chunks={{'time': {args.chunk_time}}})")
