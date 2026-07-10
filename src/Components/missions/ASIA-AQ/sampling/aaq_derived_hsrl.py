#!/usr/bin/env python3
"""
    Derives parameters that can be compared to AAQ HSRL observations
"""
import os, sys
import xarray as xr
import argparse
import yaml
from glob import glob
import numpy as np
from pyobs.aop import G2GAOP
from pyobs.sampler import addVertCoord

def calc_hsrl_pm25(traj_ds,config,h5f):

    # ----------------
    # calculate HSRL PM2.5
    # ----------------
    optics = G2GAOP(traj_ds,config=config['hsrl_pm25_config'])

    # total PM at 35% RH
    pm = optics.getPM(pmsize=2.5,aerodynamic=True,fixrh=0.35)

    # speciated PM
    for spc in optics.mieTable:
        pm_spc = optics.getPM(pmsize=2.5,Species=spc,aerodynamic=True,fixrh=0.35)
        pm = pm.assign({spc:pm_spc['PM']})
        pm = pm.assign({spc+'_FWATER':pm_spc['FWATER']})

    # write HSRL PM25 to netcdf file
    pm_hsrl = xr.Dataset(pm)
    pm_hsrl['Z'] = traj_ds['Z']
    pm_hsrl['DZ'] = traj_ds['DZ']
    pm_hsrl['PS'] = traj_ds['PS']
    outFile = config['sampled_outdir'] + '/' + os.path.basename(h5f)[:-3] + '_pm25_hsrl.nc4'
    pm_hsrl.to_netcdf(outFile,engine='netcdf4')

    # total PM at 0 RH
    pm = optics.getPM(pmsize=2.5,aerodynamic=True,fixrh=0.0)

    # speciated PM
    for spc in optics.mieTable:
        pm_spc = optics.getPM(pmsize=2.5,Species=spc,aerodynamic=True,fixrh=0.0)
        pm = pm.assign({spc:pm_spc['PM']})
        pm = pm.assign({spc+'_FWATER':pm_spc['FWATER']})

    # write HSRL PM25 to netcdf file
    pm_hsrl = xr.Dataset(pm)
    pm_hsrl['Z'] = traj_ds['Z']
    pm_hsrl['DZ'] = traj_ds['DZ']    
    pm_hsrl['PS'] = traj_ds['PS']
    outFile = config['sampled_outdir'] + '/' + os.path.basename(h5f)[:-3] + '_pm25_dry_hsrl.nc4'
    pm_hsrl.to_netcdf(outFile,engine='netcdf4')

def calc_hsrl_pm10(traj_ds,config,h5f):

    # ----------------
    # calculate HSRL PM10
    # ----------------
    optics = G2GAOP(traj_ds,config=config['hsrl_pm10_config'])

    # total PM at 35% RH
    pm = optics.getPM(pmsize=10.0,aerodynamic=True,fixrh=0.35)

    # speciated PM
    for spc in optics.mieTable:
        pm_spc = optics.getPM(pmsize=10.0,Species=spc,aerodynamic=True,fixrh=0.35)
        pm = pm.assign({spc:pm_spc['PM']})
        pm = pm.assign({spc+'_FWATER':pm_spc['FWATER']})

    # write HSRL PM25 to netcdf file
    pm_hsrl = xr.Dataset(pm)
    pm_hsrl['Z'] = traj_ds['Z']
    pm_hsrl['DZ'] = traj_ds['DZ']
    pm_hsrl['PS'] = traj_ds['PS']
    outFile = config['sampled_outdir'] + '/' + os.path.basename(h5f)[:-3] + '_pm10_hsrl.nc4'
    pm_hsrl.to_netcdf(outFile,engine='netcdf4')

    # total PM 0% RH
    pm = optics.getPM(pmsize=10.0,aerodynamic=True,fixrh=0.0)

    # speciated PM
    for spc in optics.mieTable:
        pm_spc = optics.getPM(pmsize=10.0,Species=spc,aerodynamic=True,fixrh=0.0)
        pm = pm.assign({spc:pm_spc['PM']})
        pm = pm.assign({spc+'_FWATER':pm_spc['FWATER']})

    # write HSRL PM25 to netcdf file
    pm_hsrl = xr.Dataset(pm)
    pm_hsrl['Z'] = traj_ds['Z']
    pm_hsrl['DZ'] = traj_ds['DZ']
    pm_hsrl['PS'] = traj_ds['PS']
    outFile = config['sampled_outdir'] + '/' + os.path.basename(h5f)[:-3] + '_pm10_dry_hsrl.nc4'
    pm_hsrl.to_netcdf(outFile,engine='netcdf4')


if __name__ == '__main__':


    parser = argparse.ArgumentParser()
    parser.add_argument("config",help='configuration yaml file')

    args = parser.parse_args()

    config = yaml.safe_load(open(args.config))

    # get g3 files
    h5Files = sorted(glob(config['g3_h5']+'/*h5'))

    for h5f in h5Files:
        print(f"+++++++ AAQ Derive on {h5f}")
        # Read in sampled Aerosol Collection
        # and do some addtional transformations
        # --------------------------------------
        ctl = config['model_aer_ctl']
        outFile = config['sampled_outdir'] + '/' + os.path.basename(h5f)[:-3] + '_' + ctl + '.nc4'
        try:
            if os.path.exists(outFile):
                traj_ds = xr.open_dataset(outFile)
                traj_ds.pipe(addVertCoord)
            else:
                print(f"Sampled file {outFile} does not exist")
                print("Skipping....")
                continue
        except:
            print(f"Problem with reading {outFile}")
            print("Skipping....")
            continue

        if config.get('do_hsrl_pm25'): calc_hsrl_pm25(traj_ds,config,h5f)
        if config.get('do_hsrl_pm10'): calc_hsrl_pm10(traj_ds,config,h5f)
