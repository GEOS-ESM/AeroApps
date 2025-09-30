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


def plot(model='c180R_arcsix.inst2d_hwl_x', varn='TOTEXTTAU550',station='Pituffik',
         y0=0,y1=1.,scale=1.,units='',title=''):

    fp_dataset = 'stn_samples/%s.stations.2024????.nc'%(model)
    fp_dataset = 'stn_samples/%s.stations.2024.nc'%(model)
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
    for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
                 ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(20)
    plt.title('%s, %s'%(station,titlestr), size=28)
    plt.ylabel('%s %s'%(titlestr,unitstr))
    plt.savefig('stn_plots/%s.%s.%s.png'%(model,station,varn))
    plt.close(fig)

if __name__ == "__main__":
    stations = ['Ragged_Point','Capo_Verde']
    for station in stations:
        print(station)
        model   = 'c180R_arcsix.inst3d_aer_v'
        aodvar = 'TOTEXTTAU550'
        model   = 'fp.inst3_3d_aer_Nv'
        aodvar = 'TOTEXTTAU'
#        model   = 'MERRA2.tavg1_2d_aer_Nx'
#        aodvar = 'TOTEXTTAU'
        plot(model=model,station=station,varn=aodvar,units='550 nm',title='AOD')
        plot(model=model,station=station,varn='DUSMASS',scale=1.e9,units='ug m-3',title='Dust Surface Mass Concentration')
#        plot(model=model,station=station,varn='PM',scale=1.e9,units='ug m-3',title='Particulate Matter')
#        plot(model=model,station=station,varn='PLS',scale=86400.,units='mm day-1',title='Precipitation')
#        plot(model=model,station=station,varn='SLP',scale=0.01,units='hPa',title='Sea Level Pressure')

