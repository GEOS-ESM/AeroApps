#!/usr/bin/env python3
"""
Create a virtual dataset for a GEOS collection
"""

import fsspec
import json
from fsspec.implementations.reference import LazyReferenceMapper
import os, sys
from pyobs.xrctl import parse_ctl
import yaml
import argparse
from datetime import datetime
import kerchunk.hdf
import kerchunk.df
import kerchunk.combine
import dask
from dask.distributed import Client, LocalCluster

def process_single_file(path,out_name):
    """Process a single NetCDF file to create kerchunk reference"""
    try:
        with fsspec.open(path, mode='rb') as f:
            # Create references; 'path' is used as the target URL for the chunks
            scanner = kerchunk.hdf.SingleHdf5ToZarr(f, path, inline_threshold=100)
            with open(out_name, 'wb') as out_f:
                out_f.write(json.dumps(scanner.translate()).encode())
        return f"Successfully processed: {path}"
    except Exception as e:
        return f"Error processing {path}: {str(e)}"

if __name__ == '__main__':


    parser = argparse.ArgumentParser()
    parser.add_argument("config",help='configuration yaml file')
    parser.add_argument("--overwrite",action='store_true',help='overswrite existing json files')
    parser.add_argument("--append",action='store_true',help='append to existing mzz parquet')
    parser.add_argument("--n_workers",default=120,help='number of dask workers to use')

    args = parser.parse_args()

    config = yaml.safe_load(open(args.config))

    ctls = [config['model_aer_ctl']] + config['model_other_ctl']

    # create output directory
    ind_outdir = config['individual_references']
    if not os.path.exists(ind_outdir):
        os.makedirs(ind_outdir)

    com_outdir = config['combined_references']
    if not os.path.exists(com_outdir):
        os.makedirs(com_outdir)

    fs_local = fsspec.filesystem("file")

    for ctl in ctls:
        if ctl:

            # get the time_range of files on disk
            time_start = None
            time_end = None
            istart = 0
            iend = 0
            paths = parse_ctl(ctl)
            for i,p in enumerate(paths):
                if os.path.exists(p) and time_start is None:
                    time_start = datetime.strptime(os.path.basename(p).split('.')[-2],'%Y%m%d_%H%Mz')
                    time_end = time_start
                    istart= i
                    iend = i
                elif os.path.exists(p):
                    date = datetime.strptime(os.path.basename(p).split('.')[-2],'%Y%m%d_%H%Mz')
                    if date > time_end: time_end = date
                    iend = i
            
            paths = paths[istart:iend+1]
            time_range = [time_start,time_end]
            
            # Create delayed tasks for each file
            delayed_tasks = []
            json_list = []
            for path in paths:
                out_name = f"{ind_outdir}/{os.path.basename(path)}.json"
                if os.path.exists(out_name) and not args.overwrite:
                    print(f"{out_name} exists, continuing")
                    if not args.append:
                        json_list += fsspec.filesystem('file').glob(out_name)
                else:
                    delayed_tasks.append(dask.delayed(process_single_file)(path,out_name))
                    json_list += fsspec.filesystem('file').glob(out_name)

            # Execute all tasks in parallel
            if len(delayed_tasks):
                print(f"Processing {len(delayed_tasks)} files using Dask...")
                n_workers = min([args.n_workers,len(delayed_tasks)])
                kwargs = dict(n_workers=n_workers, threads_per_worker=1, memory_limit='4GB')
                with LocalCluster(**kwargs) as cluster, Client(cluster) as client:
                    results = dask.compute(*delayed_tasks)

                    # Print results
                    for result in results:
                        print(result)

            # Define the dimension you want to concatenate along (e.g., 'time')
            concat_dims = ['time'] 

            # List dimensions that are identical across files and should NOT be concatenated

            # figure out if there is a lev dimension
            f = open(ctl)
            content = f.read()
            f.close()
            if 'zdef' in content:
                identical_dims = ['lat', 'lon','lev'] 
            else:
                identical_dims = ['lat', 'lon']

#                # Collect all individual JSON paths
#                json_list = sorted(fsspec.filesystem('file').glob(f'{ind_outdir}/*.json'))

            # Combine them along the time dimension
            mzzoutparq = com_outdir + f'/mzz_{ctl}.parq'
            if args.append:
                out = LazyReferenceMapper(root=mzzoutparq, fs=fs_local)
                mzz = kerchunk.combine.MultiZarrToZarr.append(
                    path=json_list,
                    original_refs=out,
                    remote_protocol='file',
                    concat_dims=concat_dims,
                    identical_dims=identical_dims,
                    # "INDEX" tells Kerchunk to use 0, 1, 2... based on file order
                    coo_map={"time": "INDEX"}
                )
            else:
                out = LazyReferenceMapper.create(root=mzzoutparq, fs=fs_local, record_size=5000)
                mzz = kerchunk.combine.MultiZarrToZarr(
                    path=json_list,
                    remote_protocol='file',
                    concat_dims=concat_dims,
                    identical_dims=identical_dims,
                    # "INDEX" tells Kerchunk to use 0, 1, 2... based on file order
                    coo_map={"time": "INDEX"},
                    out=out
                )

            # Translate the combined references into a single dictionary object
            combined_refs = mzz.translate()

            out.flush()  #make sure everything is written out

            print(f"Successfully created combined mzz file: {mzzoutparq}")

            # Convert the dictionary to a Parquet reference store
            outparq = com_outdir + f'/combined_{ctl}_{time_start.isoformat()}_{time_end.isoformat()}.parq'
            kerchunk.df.refs_to_dataframe(
                combined_refs, 
                outparq, 
                storage_options={'compression': 'zstd'}  # Zstandard compression is highly efficient
            )


            print(f"Successfully created combined reference file: {outparq}")

