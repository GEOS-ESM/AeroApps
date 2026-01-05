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
import kerchunk.df
import kerchunk.combine
from dask.distributed import Client, LocalCluster
import dask.bag as db
import dask



if __name__ == '__main__':


    parser = argparse.ArgumentParser()
    parser.add_argument("config",help='configuration yaml file')
    parser.add_argument("--n_workers",type=int,default=120,help='number of pool workers to use')

    args = parser.parse_args()

    config = yaml.safe_load(open(args.config))

    ctls = [config['model_aer_ctl']] + config['model_other_ctl']

    # get directory of json kerchunk files
    ind_outdir = config['individual_references']

    com_outdir = config['combined_references']
    if not os.path.exists(com_outdir):
        os.makedirs(com_outdir)

    fs_local = fsspec.filesystem("file")

    for ctl in ctls:
        if ctl:
            
            # get list of json files
            json_list = sorted(fsspec.filesystem('file').glob(f'{ind_outdir}/*{ctl}*.json'))
            json_split = np.array_split(np.array(json_list),int(len(json_list)/2)) #<- create smaller lists of individual json files
            json_split = [list(x) for x in json_split]

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


            # write the first MZZ
            mzzoutparq = com_outdir + f'/mzz_{ctl}.parq'
            out = LazyReferenceMapper.create(root=mzzoutparq, fs=fs_local, record_size=500000)
            first_batch = json_split.pop(0)
            mzz = kerchunk.combine.MultiZarrToZarr(
                path=first_batch,
                remote_protocol='file',
                concat_dims=concat_dims,
                identical_dims=identical_dims,
                # "INDEX" tells Kerchunk to use 0, 1, 2... based on file order
                coo_map={"time": "INDEX"},
                out=out
            ).translate()

            out.flush() #make sure everything is written out
            time_offset = len(first_batch)

            # define a generic MZZ function
            def multi_multizarr(json_list,concat_dims,identical_dims,time_offset):
                time_coords = list(range(time_offset, time_offset + len(json_list)))
                coo_map = {"time": time_coords}
                mzz = kerchunk.combine.MultiZarrToZarr(
                    path=json_list,
                    remote_protocol='file',
                    concat_dims=concat_dims,
                    identical_dims=identical_dims,
                    coo_map=coo_map
                )


                return mzz.translate()


            # do combine the rest in batches and append
            kwargs = dict(n_workers=args.n_workers)
            with LocalCluster(**kwargs) as cluster, Client(cluster) as client:
                for istart in range(0,len(json_split),args.n_workers):

                    iend = min([istart + args.n_workers,len(json_split)])

                    current_batch = json_split[istart:iend]
        
                    # Calculate time offsets for each batch in this worker group
                    batch_info = []
                    for batch in current_batch:
                        batch_info.append((batch, concat_dims, identical_dims, time_offset))
                        time_offset += len(batch)


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
                        # "INDEX" tells Kerchunk to use 0, 1, 2... based on file order
#                        coo_map={"time": "INDEX"}
                    ).translate()

                    out.flush()  #make sure everything is written out
                    client.wait_for_workers(args.n_workers)

            print(f"Successfully created combined mzz file: {mzzoutparq}")

            # Read the complete accumulated references from the LazyReferenceMapper
            final_out = LazyReferenceMapper(root=mzzoutparq, fs=fs_local)

            # Convert the dictionary to a Parquet reference store
            outparq = com_outdir + f'/reference_store_{ctl}.parq'
            kerchunk.df.refs_to_dataframe(
                final_out, 
                outparq, 
                storage_options={'compression': 'zstd'}  # Zstandard compression is highly efficient
            )


            print(f"Successfully created combined reference file: {outparq}")

