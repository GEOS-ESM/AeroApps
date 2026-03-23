#!/usr/bin/env python3
"""
    Script to create plots comparing a GEOS simulations to each other sampled at IMPROVE sites 
    Reads in the outputs from sample/improve_sampler.py
"""
import os, sys
import argparse
import yaml
from pyobs.improve import SITE_MAP, IMPROVE
import xarray as xr
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats
from ts_plot_gas import ts_plot

def clean_node(ds):
    if ds is None:
        return ds
    # Add .copy() to get a mutable version
    ds = ds.copy()
    
    for var in list(ds.variables):
        if ds[var].dtype == object:
            print(f'  Converting {var} to string')
            try:
                # Try converting to string
                ds[var] = ds[var].astype(str)
            except Exception as e:
                # If conversion fails, drop the variable
                print(f'  Could not convert {var}, dropping: {e}')
                ds = ds.drop_vars(var)
    
    return ds

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("config",help='configuration yaml file')

    args = parser.parse_args()

    config = yaml.safe_load(open(args.config))

    # get model sampled data 
    outdir = config['sampled_outdir']
    baseline_ctl    = config['model_baseline']
    oxh_ctl         = config['model_oxh']
    oxm_ctl         = config['model_oxm']
    spcs = ['NH3','HNO3CONC']
    units = ['kg kg-1','kg m-3']

    # get output directory
    plot_outdir = config['plot_outdir'] + '/' + baseline_ctl
    os.makedirs(plot_outdir,exist_ok=True)

    baseline_ds = xr.open_dataset(f'{outdir}/improve.{baseline_ctl}.nc4')
    oxh_ds      = xr.open_dataset(f'{outdir}/improve.{oxh_ctl}.nc4')
    oxm_ds      = xr.open_dataset(f'{outdir}/improve.{oxm_ctl}.nc4')

    for spc,unit in zip(spcs,units):
        # Get baseline data
        ds_baseline = baseline_ds[spc]
        ds_oxh = oxh_ds[spc]
        ds_oxm = oxm_ds[spc]
        
        
        # Create dataframes for baseline vs obs (top subplot)
        monthly_mean = ds_baseline.groupby('time.month').mean(skipna=True).rename({'month': 'time'})
        df_baseline = monthly_mean.to_dataframe(name=spc).reset_index()
        df_baseline["Source"] = "Baseline"
        
        monthly_mean = ds_oxh.groupby('time.month').mean(skipna=True).rename({'month': 'time'})
        df_oxh = monthly_mean.to_dataframe(name=spc).reset_index()
        df_oxh["Source"] = 'OxH'

        monthly_mean = ds_oxm.groupby('time.month').mean(skipna=True).rename({'month': 'time'})
        df_oxm = monthly_mean.to_dataframe(name=spc).reset_index()
        df_oxm["Source"] = 'OxM'
        
        # Create dataframe for top subplot (baseline vs obs)
        plot_df_top = pd.concat([df_baseline, df_oxh, df_oxm], ignore_index=True)
        
        
        # Create the plot
        ts_plot(plot_df_top, spc, unit, plot_outdir)

         
