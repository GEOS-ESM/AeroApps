#!/usr/bin/env python3
"""
Create a virtual dataset for a GEOS collection
"""

import fsspec
from fsspec.implementations.reference import LazyReferenceMapper
import json
import os, sys
from pyobs.xrctl import parse_ctl
import yaml
import argparse
import numpy as np
from datetime import datetime
import kerchunk.hdf
import kerchunk.df
import kerchunk.combine
from dask.distributed import Client, LocalCluster
import dask.bag as db
import dask
import xarray as xr
import multiprocessing

if __name__ == '__main__':


    parser = argparse.ArgumentParser()
    parser.add_argument("config",help='configuration yaml file')
    parser.add_argument("--append",action='store_true',help='append to an existing parquet file')
    parser.add_argument("--n_workers",type=int,default=120,
                        help='number of pool workers to use (default=120)')

    args = parser.parse_args()

    config = yaml.safe_load(open(args.config))

    ctls = config['reference_ctl_list']
    if len(ctls) == 0:
        raise ValueError('No reference control files given in configuration file. Cannot proceed.')

    com_outdir = config['combined_references_dir']
    os.makedirs(com_outdir,exist_ok=True)

    fs_local = fsspec.filesystem("file")

    for ctl in ctls:
        # get the list files on disk
        paths = parse_ctl(ctl)
        path_list = []
        for p in paths:
            if os.path.exists(p): path_list.append(p)
        path_list = np.array(path_list)

        if len(path_list) == 0:
            raise ValueError(f'No {ctl} files found on disk. Cannot proceed.')

        # get the name of the reference store
        mzzoutparq = com_outdir + f'/reference_store_{ctl}.parq'

        # if appending, figure out which files to add to existing reference store
        if args.append:
            # figure out last date in parquet
            fs = fsspec.filesystem("reference", fo=mzzoutparq, remote_protocol='file', lazy=True)
            ds = xr.open_dataset(fs.get_mapper(), engine="zarr",chunks='auto', backend_kwargs={"consolidated": False})
            last_time = ds.time[-1].values
            # convert to datetime
            last_time = last_time.astype('M8[us]').astype('O')


            # figure out which files come after this last time
            date_list = []
            for p in path_list:
                date_list.append(datetime.strptime(os.path.basename(p).split('.')[-2],'%Y%m%d_%H%Mz'))

            date_list = np.array(date_list)
            keep = date_list > last_time
            path_list = path_list[keep]
            if len(path_list) == 0:
                raise ValueError(f'No {ctl} files found to append on disk. Cannot proceed.')


        # define function to convert HDF5 to Zarr
        def process_single_file(path):
            """Process a single NetCDF file to create kerchunk reference"""
            try:
                with fsspec.open(path, mode='rb') as f:
                    # Create references; 'path' is used as the target URL for the chunks
                    scanner = kerchunk.hdf.SingleHdf5ToZarr(f, path, inline_threshold=0).translate()
                return scanner
            except Exception as e:
                raise RuntimeError(f"Error processing {path}: {str(e)}") 

 
        # split the list of files into chunks
        path_split = np.array_split(path_list,int(len(path_list)/2)) #<- create smaller lists of individual files
        path_split = [list(x) for x in path_split]

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

        if not args.append:
            # write the first MZZ
            out = LazyReferenceMapper.create(root=mzzoutparq, fs=fs_local, record_size=100000)
            first_batch = path_split.pop(0)

            n_workers = min([args.n_workers,len(first_batch)])
            with multiprocessing.Pool(n_workers) as pool:
                single_refs = pool.map(process_single_file, first_batch)            

            # 'cf:time' decodes internal 'time' variable using CF-conventions
            coo_map={"time": "cf:time"}
            # 'M8[ns]' forces the result to NumPy datetime64 nanoseconds
            coo_dtypes={"time": "M8[ns]"}

            mzz = kerchunk.combine.MultiZarrToZarr(
                path=single_refs,
                remote_protocol='file',
                concat_dims=concat_dims,
                identical_dims=identical_dims,
                coo_map=coo_map,
                coo_dtypes=coo_dtypes,
                out=out
            ).translate()

            out.flush() #make sure everything is written out

        # define a generic MZZ function
        def multi_multizarr(single_refs,concat_dims,identical_dims):
            coo_map = {"time": "cf:time"}
            coo_dtypes={"time": "M8[ns]"}
            mzz = kerchunk.combine.MultiZarrToZarr(
                path=single_refs,
                remote_protocol='file',
                concat_dims=concat_dims,
                identical_dims=identical_dims,
                coo_map=coo_map,
                coo_dtypes=coo_dtypes
            )


            return mzz.translate()


        # combine the rest in batches and append
        kwargs = dict(n_workers=args.n_workers)
        with LocalCluster(**kwargs) as cluster, Client(cluster) as client:
            for istart in range(0,len(path_split),args.n_workers):

                iend = min([istart + args.n_workers,len(path_split)])

                current_batch = path_split[istart:iend]
                flat_batch = [item for sublist in current_batch for item in sublist]

                intermediates = db.from_sequence(flat_batch)
                single_refs   = intermediates.map(process_single_file).compute()

                singles_batch = np.array_split(np.array(single_refs),len(current_batch))
                singles_batch = [list(x) for x in singles_batch]
                # upload the data to the workers' memory once and return Dask Futures
                singles_batch = [client.scatter(batch) for batch in singles_batch]

                # package together inputs for each batch
                batch_info = []
                for batch in singles_batch:
                    batch_info.append((batch, concat_dims, identical_dims))


                intermediates = db.from_sequence(batch_info)
                results = intermediates.starmap(multi_multizarr).compute()

                # append to existing out
                out = LazyReferenceMapper(root=mzzoutparq, fs=fs_local)
                mzz = kerchunk.combine.MultiZarrToZarr.append(
                    path=results,
                    original_refs=out,
                    remote_protocol='file',
                    concat_dims=concat_dims,
                    identical_dims=identical_dims,
                ).translate()

                out.flush()  #make sure everything is written out
                client.wait_for_workers(args.n_workers)

        print(f"Successfully created combined mzz file: {mzzoutparq}")



