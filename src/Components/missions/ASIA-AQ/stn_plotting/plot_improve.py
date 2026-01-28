#!/usr/bin/env python3
"""
    Script to create plots comparing a GEOS simulation to IMPROVE observations
    Reads in the outputs from sample/improve_derived.py
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

def ts_plot(plot_df,spcname):

    plt.figure(figsize=(12, 6))

    # Plot both lines with 25th-75th percentile shading
    # ("pi", 50) creates an interval covering 50% of the data (25th to 75th)
    sns.lineplot(
        data=plot_df,
        x="time",
        y="PM25",
        hue="Source",
        errorbar=("pi", 50), 
        alpha=0.8
    )

    plt.title(f"Comparison of Surface {spcname} PM2.5: GEOS vs. IMPROVE (Mean & 25-75th Percentile)")
    plt.ylabel(r"Concentration ($\mu g/m^3$)")
    plt.xlabel("Date")
    plt.grid(True, alpha=0.3)
    plt.show()

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("config",help='configuration yaml file')

    args = parser.parse_args()

    config = yaml.safe_load(open(args.config))

    # get model sampled data 
    outdir = config['sampled_outdir']
    ctl    = config['model_aer_ctl']
    spcs = ['BC','BR','DU','NH4','NI','OC','SS','SU','TOTAL']
    stn_ds = {}
    for spc in spcs:
        stn_ds[spc] = xr.open_dataset(f'{outdir}/improve.{ctl}.pm25_{spc}.nc4')

    # get time span of data
    syear = pd.Timestamp(stn_ds['TOTAL'].time[0].values).year
    eyear = pd.Timestamp(stn_ds['TOTAL'].time[-1].values).year

    yearlist = list(range(syear,eyear+1))

    # get IMPROVE site locations
    site_path = config['improve_site_map']
    sites = SITE_MAP(site_path)

    # get IMPROVE data
    imp_path = config['improve_data_dir']
    implist = []
    for year in yearlist:
        implist += [f'{imp_path}/IMPAER_{year}.txt']
    imp = IMPROVE(implist,site_path=site_path)

    # convert IMPROVE dataframe to xarray dataset for easier comparisons
    df = imp.df.rename(columns={'SiteCode':'station','Start_Time':'time','Longitude':'lon','Latitude':'lat'})
    df = df.set_index(['time','station','ParamCode'])
    print(f"Original length: {len(df.index)}")
    print(f"Unique index length: {len(df.index.unique())}")
    # there  can be duplicate entries because there are multiple instruments at one site
    # let's just keep the first measurement in the database
    df_unique = df[~df.index.duplicated(keep='first')]
    # Now convert the cleaned DataFrame to xarray
    imp_ds = df_unique.to_xarray()
    # set all -999 to nan
    imp_ds = imp_ds.where(imp_ds != -999.)
    # retain all pm25 values to where Status == 'V0', set all the rest to nan
    imp_ds["pm25"] = imp_ds["pm25"].where(imp_ds["Status"] == "V0")


    # subset geos and improve for the dates that are available in both
    for spc in spcs:
        geos_sync, imp_sync = xr.align(stn_ds[spc], imp_ds, join="inner")
        stn_ds[spc] = geos_sync

    stn_ds['improve'] = imp_sync
   
    # create a datatree of the two datasets
    dt = xr.DataTree.from_dict(stn_ds) 

    ## TOTAL PM2.5
    imspc ={'RCFM': dt['TOTAL'], 
            'SeaSaltf': dt['SS'], 
            'SOILf': dt['DU'],
            'ECf': dt['BC'],
            'OMCf': dt['OC'] + dt['BR'],
            'ammSO4f': dt['SU']
            }
    for spc in imspc:
        df_geos = imspc[spc].PMDRY.to_dataframe(name="PM25").reset_index()
        df_geos["Source"] = "GEOS"
        df_imp   = dt['improve'].pm25.sel(ParamCode=spc,drop=True).to_dataframe(name="PM25").reset_index()
        df_imp["Source"] = 'IMPROVE'
        plot_df = pd.concat([df_geos, df_imp],ignore_index=True) 
   
        ts_plot(plot_df,spc)
     
