#!/usr/bin/env python3

import pandas as pd
import xarray as xr
from pyobs.sampler import STATION
from glob import glob
from datetime import datetime,timedelta
import sys, os
import matplotlib.pyplot as plt
import numpy as np
import argparse
if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("station_file",
           help="input CSV file with starting points name,lon,lat")

    args = parser.parse_args()

    # read start point CSV
    # ---------------------
    station = pd.read_csv(args.station_file)
    stations = station.name.values[0]
   
    ctlFile = 'fp/opendap/seamless/tavg3_2d_aer_Nx.latest'
    plt_Variables= ['TOTEXTTAU','DUEXTTAU','SSEXTTAU','SUEXTTAU','PM25','CCEXTTAU']

    # start of timeseries
    sd = datetime(2024,1,1,1,30)

    # end of time series
    # today's date + 5 day forecast, rounded to closest synoptic mid time
    today = datetime.today()
    ed = today + timedelta(days=4)
    hour = 3*(ed.hour//3)
    ed = datetime(ed.year,ed.month,ed.day,ed.hour+1,30)

    # create a directory for figures
    opath = 'plots_{}'.format(today.strftime('%Y%m%d'))
    if not os.path.exists(opath):
        os.makedirs(opath)

    # read in climatology
    stats = ['mean','med','min','max','25','75']
    clim = {}
    for stat in stats:
        dataset = sorted(glob('{}_clim_{}.nc'.format(stations,stat)))
        clim[stat] = xr.open_mfdataset(dataset,parallel=True)

    Variables = list(clim['mean'].keys())
    Variables.remove('PM25')
    Variables.remove('CCEXTTAU')

    # read seamless control files
    f = open(ctlFile,'r')
    dsetRoot = f.readline().split()[1][:-4]
    ch, assim_s, assim_e, assim_path = f.readline().split()
    ch, forec_s, forec_e, forec_path = f.readline().split()
    for l in range(7):
        f.readline()
    tdef,t_e,linear,ctl_start, dmin = f.readline().split()
    ctlsdate = datetime.strptime(ctl_start,'%H:%MZ%d%b%Y')
    
    assim_sdate = datetime.strptime(ctl_start,'%H:%MZ%d%b%Y')
    assim_edate = assim_sdate + timedelta(hours=(3*int(assim_e)))

    forec_sdate = assim_sdate + timedelta(hours=(3*int(forec_s)))

    f.close()

    # get FP filenames
    fp_dataset = []
    dt = timedelta(hours=3)
    tsd = sd
    while sd <= ed:
        yy = sd.strftime('%Y')
        mm = sd.strftime('%m')
        dd = sd.strftime('%d')
        hh = sd.strftime('%H')
        nn = sd.strftime('%M')
        if sd <= assim_edate:
            filename = assim_path.replace('%y4',yy).replace('%m2',mm).replace('%d2',dd).replace('%h2',hh).replace('%n2',nn)
        else:
            filename = forec_path.replace('%y4',yy).replace('%m2',mm).replace('%d2',dd).replace('%h2',hh).replace('%n2',nn)
        fp_dataset.append(dsetRoot + filename)

        sd += dt

    # open dataset
    # need the preprocess because forecast has additional variables
    # without it open_mfdataset takes FOREVER
    ds = xr.open_mfdataset(fp_dataset,parallel=True,preprocess=lambda ds: ds[Variables])
    ds['PM25'] = ds['BCSMASS']+ds['OCSMASS']+ds['DUSMASS25']+ds['SO4SMASS']+ds['SSSMASS25']
    ds['CCEXTTAU'] = ds['OCEXTTAU']+ds['BCEXTTAU']

    # initialize stn and sample
    stn = STATION(station.name,station.lon,station.lat,ds,verbose=True)
    ds = stn.sample()
    
    # start plotting stuff
    for var in plt_Variables:
        # load so plotting is faster
        ds[var].load()

        fig = plt.figure(figsize=(13,6))
        ax = fig.add_subplot(111)

        ax.fill_between(clim['min'].time.values,
            clim['min'][var].values.squeeze(),
            clim['max'][var].values.squeeze(),
            alpha=0.2,
            label='min-max climatology')

        ax.fill_between(clim['25'].time.values,
            clim['25'][var].values.squeeze(),
            clim['75'][var].values.squeeze(),
            alpha=0.2,
            label='q25-q75 climatology')

        ax.plot(clim['med'].time.values,clim['med'][var].values.squeeze(),'k-',
            label='median climatology')

        iassim = ds.time.values <= np.datetime64(assim_edate)
        ax.plot(ds.time.values[iassim],ds[var].values.squeeze()[iassim],'b-',label='2024 assimilation')
        ifor = ds.time.values >= np.datetime64(assim_edate)
        ax.plot(ds.time.values[ifor],ds[var].values.squeeze()[ifor],'r-',label='2024 forecast')

#        ylim = ax.get_ylim()
#        ax.plot([today,today],ylim,'g-',label='today')

        ax.set_title('{} {}'.format(stations.upper(),var))

        if 'TAU' in var:
            title = 'Aerosol Optical Depth'
        else:
            title = 'Aerosol Surface Mass [kg m-3]'

        ax.set_ylabel(title)
        ax.tick_params(axis='x', labelrotation=25)

        plt.legend()
        plt.savefig('{}/{}_{}_{}.png'.format(opath,stations,var,today.strftime('%Y%m%d_%H')))
#        plt.show()
        plt.close(fig)
    





