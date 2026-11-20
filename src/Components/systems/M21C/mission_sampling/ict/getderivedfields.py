#!/usr/bin/env python3
"""
    Derives parameters that can be compared to airborne observations
"""
import os, sys
import xarray as xr
import argparse
import yaml
from glob import glob
import numpy as np
from pyobs.aop import G2GAOP

def calc_ams(traj_ds,config):
    print(config['ams_config'])
    # ----------------
    # calculate AMS PM
    # ----------------
    optics = G2GAOP(traj_ds,config=config['ams_config'])

    # total PM
    pm = optics.getPM(pmsize=config['pmsize'],vacuum_aerodynamic=True,fixrh=0.0)

    # speciated PM
    for spc in optics.mieTable:
        pm_spc = optics.getPM(pmsize=config['pmsize'],Species=spc,vacuum_aerodynamic=True,fixrh=0.0)
        pm = pm.assign({spc:pm_spc['PM']})


    # AMS observes at STP (273 K & 1013 mb)
    # so need to convert ambient concentration to STP
    # using the density of Air at STP = 1.2754 kg/m3
    pm_ams = {}
    for spc in ['PM'] + list(optics.mieTable.keys()):
        pm_stp  = pm[spc]*(1.2754/pm['AIRDENS'])
        pm_stp.attrs.update(pm[spc].attrs)
        pm_ams[spc] = pm_stp

    # MSA doesn't have an optics file, so need to do it MSA seperately.
    # Assume total mass is below 1.5 micron
    q_conc = (traj_ds['AIRDENS'] * traj_ds['MSA']*1e9)  # Aerosol mass concentration in ug/m3
    pm_ams['MSA'] = q_conc*(1.2754/pm['AIRDENS'])
    attrs = {'long_name':'Particulate Matter', 'units':'microgram m-3'}
    pm_ams['MSA'].attrs.update(attrs)

    # write AMS PM to netcdf file
    pm_ams = xr.Dataset(pm_ams)
    if 'H' in traj_ds: #MERRA-2 does not have H in aer_Nv files
        pm_ams['H'] = traj_ds['H']  
    pm_ams['PS'] = traj_ds['PS']
    pm_ams['delp'] = pm['DELP']
    pm_ams['AIRDENS'] = pm['AIRDENS']
    outFile = config['sampled_outdir'] + '/' + os.path.basename(ict)[:-3] + 'pm_ams.nc4'
    pm_ams.to_netcdf(outFile,engine='netcdf4')

def calc_sp2(traj_ds,config):
    # --------------------
    # calculate SP2 BC Mass
    # --------------------
    optics = G2GAOP(traj_ds,config=config['sp2_config'])
    pm = optics.getPM(pmsize=config['sp2size'],fixrh=0.0)

    # SP2 observes at STP (273 K & 1013 mb)
    # so need to convert ambient concentration to STP
    # using the density of Air at STP = 1.2754 kg/m3
    pm_sp2 = {}
    for spc in ['PM']:
        pm_stp = pm[spc]*(1.2754/pm['AIRDENS'])
        pm_stp.attrs.update(pm[spc].attrs)
        pm_sp2['BC'] = pm_stp

    # write SP2 PM to netcdf file
    pm_sp2 = xr.Dataset(pm_sp2)
    if 'H' in traj_ds: #MERRA-2 does not have H in aer_Nv files
        pm_sp2['H'] = traj_ds['H']  
    pm_sp2['PS'] = traj_ds['PS']
    pm_sp2['delp'] = pm['DELP']
    pm_sp2['AIRDENS'] = pm['AIRDENS']
    outFile = config['sampled_outdir'] + '/' + os.path.basename(ict)[:-3] + 'bc_sp2.nc4'
    pm_sp2.to_netcdf(outFile,engine='netcdf4')

def calc_large_aop(traj_ds,config):
    # --------------------
    # calculate LARGE AOPs
    # --------------------

    # LARGE reports at STP
    # AIRDENS at STP = 1.2754 kg/m3
    traj_STP = traj_ds.copy()
    traj_STP['AIRDENS'][:] = 1.2754


    # BULK AEROSOLS
    # ----------------
    optics = G2GAOP(traj_STP,config=config['large_config'])
    wavs = ['450','470','532','550','660','700']
    varlist = ['EXT','SCA','BSC','DEPOL']

    # loop through wavelengths
    for wl in wavs:
        # model RH
        aop = optics.getAOPext(wavelength=float(wl),doaback=False)

        for spc in optics.mieTable:
            aop_spc = optics.getAOPext(wavelength=float(wl),Species=spc,doaback=False)
            for var in varlist:
                aop = aop.assign({var+'_'+spc:aop_spc[var]})

        aop['H'] = traj_STP['H']
        aop['PS'] = traj_STP['PS']
        outFile = config['sampled_outdir'] + '/' + os.path.basename(ict)[:-3] + f'{wl}_aop_STP.nc4'
        aop.to_netcdf(outFile,engine='netcdf4')

        # loop through RH = 0,20,40,80
        for rh in [0.0,0.2,0.4,0.8]:    
            aop = optics.getAOPext(wavelength=float(wl),doaback=False,fixrh=rh)

            for spc in optics.mieTable:
                aop_spc = optics.getAOPext(wavelength=float(wl),Species=spc,doaback=False,fixrh=rh)
                for var in varlist:
                    aop = aop.assign({var+'_'+spc:aop_spc[var]})

            aop['H'] = traj_STP['H']
            aop['PS'] = traj_STP['PS']
            rhi = int(rh*10)
            outFile = config['sampled_outdir'] + '/' + os.path.basename(ict)[:-3] + f'{wl}_aop_STP_rh{rhi}.nc4'
            aop.to_netcdf(outFile,engine='netcdf4')

def calc_large_submicron_aop(traj_ds,config):

    # calculate LARGE SUBMICRON AOPs
    # ------------------------------

    # LARGE reports at STP
    # AIRDENS at STP = 1.2754 kg/m3
    traj_STP = traj_ds.copy()
    traj_STP['AIRDENS'][:] = 1.2754

    optics = G2GAOP(traj_STP,config=config['large_submicron_config'])
    wavs = ['450','470','532','550','660','700']
    varlist = ['EXT','SCA','BSC','DEPOL']

    # get submicron aerosol mass
    pm = optics.getPM(pmsize=1.0,aerodynamic=True,fixrh=0.0)

    for spc in list(optics.mieTable.keys()):
        pm_spc = optics.getPM(pmsize=1.0,Species=spc,aerodynamic=True,fixrh=0.0)

        pm = pm.assign({spc:pm_spc['PM']*1e-9/pm_spc['AIRDENS']})  # back to kg/kg

    pm['RH'] = traj_STP['RH'].copy()
    optics.aer = pm

    # loop through wavelengths
    for wl in wavs:

        # model RH
        aop = optics.getAOPext(wavelength=float(wl),doaback=False)

        for spc in optics.mieTable:
            aop_spc = optics.getAOPext(wavelength=float(wl),Species=spc,doaback=False)
            for var in varlist:
                aop = aop.assign({var+'_'+spc:aop_spc[var]})

        aop['H'] = traj_STP['H']
        aop['PS'] = traj_STP['PS']
        outFile = config['sampled_outdir'] + '/' + os.path.basename(ict)[:-3] + f'{wl}_submicron_aop_STP.nc4'
        aop.to_netcdf(outFile,engine='netcdf4')        

        # loop through RH = 0,20,40,80
        for rh in [0.0,0.2,0.4,0.8]:

            aop = optics.getAOPext(wavelength=float(wl),doaback=False,fixrh=rh)

            for spc in optics.mieTable:
                aop_spc = optics.getAOPext(wavelength=float(wl),Species=spc,doaback=False,fixrh=rh)
                for var in varlist:
                    aop = aop.assign({var+'_'+spc:aop_spc[var]})

            aop['H'] = traj_STP['H']
            aop['PS'] = traj_STP['PS']
            irh = int(10*rh)
            outFile = config['sampled_outdir'] + '/' + os.path.basename(ict)[:-3] + f'{wl}_submicron_aop_STP_rh{irh}.nc4'
            aop.to_netcdf(outFile,engine='netcdf4')


if __name__ == '__main__':


    parser = argparse.ArgumentParser()
    parser.add_argument("config",help='configuration yaml file')

    args = parser.parse_args()

    config = yaml.safe_load(open(args.config))

    # get aircraft files
    ictFiles = sorted(glob(config['mergefiles']+'/*ict'))
    if 'flight_dates' in config and config['flight_dates']:
        valid_dates = config['flight_dates']
        ictFiles = [ict for ict in ictFiles if any(date in os.path.basename(ict) for date in valid_dates)]

    for ict in ictFiles:
        print(f"+++++++ AAQ Derive on {ict}")
        # Read in sampled Aerosol Collection
        # and do some addtional transformations
        # --------------------------------------
        ctl = config['model_aer_ctl']
        outFile = config['sampled_outdir'] + '/' + os.path.basename(ict)[:-3] + ctl + '.nc4'
        try:
            if os.path.exists(outFile):
                traj_ds = xr.open_dataset(outFile)
            else:
                print(f"Sampled file {outFile} does not exist")
                print("Skipping....")
                continue
        except:
            print(f"Problem with reading {outFile}")
            print("Skipping....")
            continue

        if config['do_ams']: calc_ams(traj_ds,config)
        if config['do_sp2']: calc_sp2(traj_ds,config)
        if config['do_large_aop']: calc_large_aop(traj_ds,config) 
        if config['do_large_submicron_aop']: calc_large_submicron_aop(traj_ds,config)

