#!/usr/bin/env python3
"""
    Sampler for ASIA-AQ obs and modeling experiments
"""
import os, sys
import subprocess
import shutil
import xarray as xr
import argparse
import yaml
from glob import glob
import numpy as np
from pyobs.sampler import TRAJECTORY
from pyobs.icartt import ICARTT

if __name__ == '__main__':


    parser = argparse.ArgumentParser()
    parser.add_argument("config",help='configuration yaml file')

    config = yaml.safe_load(open(args.config))

    # get dc-8 files
    ictFiles = sorted(glob(config['dc8_merge']+'/*'))

    sys.exit()

    # loop through icartt files
    for ict in ictFiles:
        m = ICARTT(ictFile)
        lon, lat, tyme = m.Nav['Longitude'], m.Nav['Latitude'], m.Nav['Time']

        for ctl in config['model_ctl']:
            traj = TRAJECTORY(tyme,lon,lat,ctl)

            traj_ds = traj.sample()
            traj_ds.to_netcdf(outFile)
