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
from pyobs.aop import G2GAOP

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
    for ict in ictFiles[0:1]:
        m = ICARTT(ict)
        lon, lat, tyme = m.Longitude_BENNETT, m.Latitude_BENNETT, m.Nav['Time']

        # Sample Aerosol Collection
        # --------------------------
        ctl = config['model_aer_ctl'] 
        chunks = {'time':1, 'lev':-1, 'lat':-1, 'lon': -1}
        traj = TRAJECTORY(tyme,lon,lat,ctl,verbose=True,chunks=chunks)
#        traj_ds = traj.sample()
#        traj_ds = traj_ds.compute()

        # write out the native sampled model fields
        outFile = config['sampled_outdir'] + '/' + os.path.basename(ict)[:-3] + ctl + '.nc4'
#        traj_ds.to_netcdf(outFile,engine='netcdf4')
        traj_ds = xr.open_dataset(outFile)

        # calculate AMS PM
        # ----------------
        optics = G2GAOP(traj_ds,config=config['ams_config'])

        # total PM
        pm = optics.getPM(pmsize=1.5,vacuum_aerodynamic=True,fixrh=0.0)

        # speciated PM
        for spc in optics.mieTable:
            pm_spc = optics.getPM(pmsize=1.5,Species=spc,vacuum_aerodynamic=True,fixrh=0.0)
            pm = pm.assign({spc:pm_spc['PM']})
            

        # AMS observes at STP (273 K & 1013 mb)
        # so need to convert ambient concentration to STP 
        # using the density of Air at STP = 1.2754 kg/m3
        pm_ams = {}
        for spc in ['PM'] + list(optics.mieTable.keys()):
            pm_stp  = pm[spc]*(1.2754/pm['AIRDENS'])
            pm_stp.attrs.update(pm[spc].attrs)
            pm_ams[spc] = pm_stp

        # Do MSA seperately.  Assume total  is below 1.5 micron
        q_conc = (traj_ds['AIRDENS'] * traj_ds['MSA']*1e9)  # Aerosol mass concentration in ug/m3
        pm_ams['MSA'] = q_conc*(1.2754/pm['AIRDENS'])
        attrs = {'long_name':'Particulate Matter', 'units':'microgram m-3'}
        pm_ams['MSA'].attrs.update(attrs)

        # write AMS PM to netcdf file
        pm_ams = xr.Dataset(pm_ams)
        pm_ams['H'] = traj_ds['H']
        pm_ams['PS'] = traj_ds['PS']
        pm_ams['delp'] = pm['DELP']
        pm_ams['AIRDENS'] = pm['AIRDENS']
        outFile = config['sampled_outdir'] + '/' + os.path.basename(ict)[:-3] + 'pm_ams.nc4'
        pm_ams.to_netcdf(outFile,engine='netcdf4')
        
        # calculate AOPs
        # at 450, 470, 532, 550, 660, 700 nm
        # at model RH, RH=0, and aircraft measured RH
        # ---------------------------------------------
        optics = G2GAOP(traj_ds,config=config['large_config'])
        wavs = ['450','470','532','550','660','700']
        varlist = ['EXT','SCA','BSC','DEPOL']
        for wl in wavs:
            # wet
            aop = optics.getAOPext(wavelength=float(wl),doaback=False)

            for spc in optics.mieTable:
                aop_spc = optics.getAOPext(wavelength=float(wl),Species=spc,doaback=False)
                for var in varlist:
                    aop = aop.assign({var+'_'+spc:aop_spc[var]})

            aop['H'] = traj_ds['H']
            aop['PS'] = traj_ds['PS']
            outFile = config['sampled_outdir'] + '/' + os.path.basename(ict)[:-3] + f'{wl}_aop.nc4'
            aop.to_netcdf(outFile,engine='netcdf4')

            # dry
            aop = optics.getAOPext(wavelength=float(wl),doaback=False,fixrh=0.0)

            for spc in optics.mieTable:
                aop_spc = optics.getAOPext(wavelength=float(wl),Species=spc,doaback=False,fixrh=0.0)
                for var in varlist:
                    aop = aop.assign({var+'_'+spc:aop_spc[var]})

            aop['H'] = traj_ds['H']
            aop['PS'] = traj_ds['PS']
            outFile = config['sampled_outdir'] + '/' + os.path.basename(ict)[:-3] + f'{wl}_aop_dry.nc4'
            aop.to_netcdf(outFile,engine='netcdf4')




        

            
