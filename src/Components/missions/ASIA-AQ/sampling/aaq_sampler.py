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

if __name__ == '__main__':


    parser = argparse.ArgumentParser()
    parser.add_argument("config",help='configuration yaml file')

    args = parser.parse_args()

    config = yaml.safe_load(open(args.config))

    # get dc-8 files
    ictFiles = sorted(glob(config['dc8_merge']+'/*ict'))

    # create output directory
    if not os.path.exists(config['sampled_outdir']):
        os.makedirs(config['sampled_outdir'])
    

    # loop through icartt files
    for ict in ictFiles:
        m = ICARTT(ict)
        lon, lat, tyme = m.Longitude_BENNETT, m.Latitude_BENNETT, m.Nav['Time']

        # Sample Aerosol Collection
        # --------------------------------------
        ctl = config['model_aer_ctl'] 
        chunks = {'time':1, 'lev':-1, 'lat':-1, 'lon': -1}
        traj = TRAJECTORY(tyme,lon,lat,ctl,verbose=True,chunks=chunks)
        traj_ds = traj.sample()
        traj_ds = traj_ds.compute()

        # write out the native sampled model fields
        outFile = config['sampled_outdir'] + '/' + os.path.basename(ict)[:-3] + ctl + '.nc4'
        traj_ds.to_netcdf(outFile,engine='netcdf4')

        # If other model collections are provided
        # sample those and write out
        # -----------------------------
        for ctl in config['model_other_ctl']:
            if ctl is not None:
                traj = TRAJECTORY(tyme,lon,lat,ctl,verbose=True,chunks=chunks)
                traj_ds = traj.sample()
                traj_ds = traj_ds.compute()

                # write out the native sampled model fields
                outFile = config['sampled_outdir'] + '/' + os.path.basename(ict)[:-3] + ctl + '.nc4'
                traj_ds.to_netcdf(outFile,engine='netcdf4')


        

            
