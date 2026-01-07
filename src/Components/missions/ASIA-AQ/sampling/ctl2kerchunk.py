#!/usr/bin/env python3
"""
Create a virtual dataset for a GEOS collection
"""

import fsspec
import json
import os, sys
from pyobs.xrctl import parse_ctl
import yaml
import argparse
import kerchunk.hdf
import multiprocessing

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
    parser.add_argument("--n_workers",type=int,default=120,help='number of pool workers to use')

    args = parser.parse_args()

    config = yaml.safe_load(open(args.config))

    ctls = [config['model_aer_ctl']] + config['model_other_ctl']

    # create output directory
    ind_outdir = config['individual_references']
    if not os.path.exists(ind_outdir):
        os.makedirs(ind_outdir)

    for ctl in ctls:
        if ctl:

            # get the list files on disk
            paths = parse_ctl(ctl)
            path_list = []
            for p in paths:
                if os.path.exists(p): path_list.append(p)
            
            if path_list:
                outdir = ind_outdir + '/' + ctl
                if not os.path.exists(outdir):
                    os.makedirs(outdir)


                # Create delayed tasks for each file
                delayed_tasks = []
                json_list = []
                for path in path_list:
                    out_name = f"{outdir}/{os.path.basename(path)}.json"
                    if os.path.exists(out_name) and (os.path.getsize(out_name) > 0) and not args.overwrite:
                        print(f"{out_name} exists, continuing")
                    else:
                        delayed_tasks.append((path,out_name))
                        json_list += fsspec.filesystem('file').glob(out_name)

                # Execute all tasks in parallel
                if delayed_tasks:
                    print(f"Processing {len(delayed_tasks)} files using multiprocessing pool...")
                    n_workers = min([args.n_workers,len(delayed_tasks)])
                    with multiprocessing.Pool(n_workers) as pool:
                        results = pool.starmap(process_single_file, delayed_tasks)


                    # Print results
                    for result in results:
                        print(result)


