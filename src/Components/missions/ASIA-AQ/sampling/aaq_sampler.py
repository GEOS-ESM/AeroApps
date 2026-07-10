#!/usr/bin/env python3
"""
    Sampler for ASIA-AQ obs and modeling experiments
"""
import os, sys
import argparse
import yaml
from glob import glob
from pyobs.sampler import TRAJECTORY
from pyobs.icartt import ICARTT
from pyobs.xrctl import parse_ctl
from dask.distributed import performance_report
from dask.distributed import Client, LocalCluster
if __name__ == '__main__':


    parser = argparse.ArgumentParser()
    parser.add_argument("config",help='configuration yaml file')
    parser.add_argument("--profile",action="store_true")
    parser.add_argument("--n_workers",type=int,default=120,
                        help='number of pool workers to use (default=120)')
    parser.add_argument("--memory_limit",default='4GB',
                        help='memory limit per worker (default=4GB)')

    args = parser.parse_args()

    config = yaml.safe_load(open(args.config))

    # create output directory
    if not os.path.exists(config['sampled_outdir']):
        os.makedirs(config['sampled_outdir'])

    # get dc-8 files
    ictFiles = sorted(glob(config['dc8_merge']+'/*ict'))

    kwargs = dict(n_workers=args.n_workers, threads_per_worker=1, memory_limit=args.memory_limit)
    with LocalCluster(**kwargs) as cluster, Client(cluster) as cluster:
        client.wait_for_workers(args.n_workers)
        for ict in ictFiles:
            print(f"+++++++ AAQ Sampling on {ict}")
            m = ICARTT(ict)
            lon, lat, tyme = m.Longitude_BENNETT, m.Latitude_BENNETT, m.Nav['Time']

            # Sample Aerosol Collection
            # --------------------------------------
            ctl = config['model_aer_ctl'] 
            cs  = config.get("cs",False)

            chunks = {'time':1, 'lev':-1, 'lat':-1, 'lon': -1}
            # check that all files exist
            time_range = [tyme.min(),tyme.max()]
            paths = parse_ctl(ctl,time_range)
            exists = True
            for p in paths:
                if not os.path.exists(p): exists = False
            
            if exists:
                traj = TRAJECTORY(tyme,lon,lat,ctl,verbose=True,chunks=chunks,engine='h5netcdf')
                
                traj_ds = traj.sample()

                outFile = config['sampled_outdir'] + '/' + os.path.basename(ict)[:-3] + ctl + '.nc4'
                # Check if the file already exists, and clobber (delete) it if it does
                if os.path.exists(outFile):
                    os.remove(outFile)

                if args.profile:
                    with performance_report(filename="dask-report.html"):
                        if cs:
                            for vn in traj_ds.data_vars:
                                print(f'Computing and writing {vn}...')
                                traj_ds[vn].compute().to_netcdf(outFile, mode='a')  # append mode
                        else:
                            traj_ds = traj_ds.compute()
                            # write out the native sampled model fields
                            traj_ds.to_netcdf(outFile,engine='netcdf4')

                else:
                    if cs:
                        for vn in traj_ds.data_vars:
                            print(f'Computing and writing {vn}...')
                            traj_ds[vn].compute().to_netcdf(outFile, mode='a')  # append mode
                    else:
                        traj_ds = traj_ds.compute()
                        # write out the native sampled model fields
                        traj_ds.to_netcdf(outFile,engine='netcdf4')

                # If other model collections are provided
                # sample those and write out
                # -----------------------------
                for ctl in config['model_other_ctl']:
                    if ctl is not None:
                        try:
                            # check that all files exist
                            time_range = [tyme.min(),tyme.max()]
                            paths = parse_ctl(ctl,time_range)
                            exists = True
                            for p in paths:
                                if not os.path.exists(p): exists = False

                            if exists:
                                traj = TRAJECTORY(tyme,lon,lat,ctl,verbose=True,chunks=chunks,engine='h5netcdf')
                            else:
                                print(f"Model input files missing for {ict}")
                                print("Skipping...")
                                continue
                        except:
                            print(f"Error reading model for {ict}")
                            print("Skipping....")
                            continue

                        traj_ds = traj.sample()

                        outFile = config['sampled_outdir'] + '/' + os.path.basename(ict)[:-3] + ctl + '.nc4'
                        # Check if the file already exists, and clobber (delete) it if it does
                        if os.path.exists(outFile):
                            os.remove(outFile)

                        if cs:
                            for vn in traj_ds.data_vars:
                                print(f'Computing and writing {vn}...')
                                traj_ds[vn].compute().to_netcdf(outFile, mode='a')  # append mode
                        else:
                            traj_ds = traj_ds.compute()
                            # write out the native sampled model fields
                            traj_ds.to_netcdf(outFile,engine='netcdf4')

            else:
                print(f"Model input files missing for {ict}")
                print("Skipping...")
                continue

        client.close()
