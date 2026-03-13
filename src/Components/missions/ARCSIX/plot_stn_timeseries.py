#!/usr/bin/env python3
import xarray as xr

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.path as mpath
from matplotlib import cm #cm is colormap
import matplotlib.dates as mdates
import matplotlib.colors as colors
import matplotlib.ticker as mticker
import matplotlib
matplotlib.use('agg')


def plot(model="fp", collection="inst3-3d-AER-Nv",
         varn="TOTEXTTAU550",station="Pituffik",
         y0=0,y1=1.,scale=1.,units="",title=""):

    modeltitle='GEOS-FP'
    modname = model
    if(model == 'fp'):
        modname = "GEOS-FP"

    dirname = f"samples/ARCSIX/sampled/stations/{modname}"
    fp_dataset = f"{dirname}/ARCSIX-{modname}-{collection}-stations_Model_2024????.nc"
    ds = xr.open_mfdataset(fp_dataset)
    stn_list = ds['station'].values.tolist()
#   Find station index
    try:
        istn = stn_list.index(station)
    except ValueError:
        print("%s not in list of: "%(station), stn_list)
    
    varv  = ds[varn]*scale
    if(varn == 'PM'):
        varn2 = 'PM2.5'
        varv2 = ds['PM25']*scale
    if(varn == 'PLS'):
        varn2 = 'PCU'
        varv2 = ds[varn2]*scale
    

    unitstr = ''
    if units != '':
        unitstr = '[%s]'%(units)

    titlestr = title
    if title == '':
        titlestr = varn
    
    fig, ax = plt.subplots(figsize=(20, 6))
    time  = ds.time.values
    ntime = ds.dims['time']

    dtFmt = mdates.DateFormatter('%m/%d') # define the formatting
    plt.gca().xaxis.set_major_formatter(dtFmt) # apply the format to the desired axis
    cf  = ax.plot(time,np.squeeze(varv[:,istn]),color='blue',linewidth=4,label=varn)
    if (varn == 'PM') or (varn == 'PLS'):
        cf  = ax.plot(time,np.squeeze(varv2[:,istn]),color='orange',linewidth=4,label=varn2)
        ax.legend()
    plt.title('%s, %s'%(station,titlestr), size=20)
    plt.ylabel('%s %s'%(titlestr,unitstr))
    plt.savefig(f"{dirname}/ARCSIX-{modname}-{collection}-stations_Model_{station}_{varn}.png")
    plt.close(fig)

if __name__ == "__main__":
    stations = ['Resolute_Bay', 'Pituffik', 'Kaffeklubben_Island', 'Eureka',
                'Alert', 'Baffin_Bay', 'Belle', 'Fram_Strait', 'Lisee',
                'Northern_Arctic_Ocean', 'Station_Nord', 'Summit_Greenland',
                'Svalbard_Zeppelin']
    for station in stations:
        print(station)
        model   = 'fp'
        aodvar = 'TOTEXTTAU'
        plot(model=model,station=station,varn=aodvar,units='550 nm',title='AOD')
#        plot(model=model,station=station,varn='PM',scale=1.e9,units='ug m-3',title='Particulate Matter')
#        plot(model=model,station=station,varn='PLS',scale=86400.,units='mm day-1',title='Precipitation')
        plot(model=model,station=station,varn='SLP',scale=0.01,units='hPa',title='Sea Level Pressure')

