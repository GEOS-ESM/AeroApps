#!/usr/bin/env python3
import xarray as xr
#from pyobs.aop import G2GAOP

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.path as mpath
from matplotlib import cm #cm is colormap
import matplotlib.dates as mdates
import matplotlib.colors as colors
import matplotlib.ticker as mticker
import matplotlib
import sys
matplotlib.use('agg')

    
config = '/home/pcolarco/silo/GMAOpyobs/src/config/m2_pm25.yaml'

#Get the HALO RGB colors
import csv
from matplotlib.colors import LinearSegmentedColormap
with open('./halo_colorbar.csv', newline='') as csvfile:
    reader = csv.DictReader(csvfile, delimiter=',',fieldnames=['r','g','b'])
    i = 0
    rgb = []
    for row in reader:
        rgb.append([int(row['r'])/255.,int(row['g'])/255.,int(row['b'])/255.])
cm = LinearSegmentedColormap.from_list(
        'my_map', rgb)


def plot(model='c180R_arcsix.inst3d_aer_v', varn='T',station='Pituffik',
         y0=250,y1=310.,scale=1.,units='',title='',species=None):

    fp_dataset = 'stn_samples/%s.stations.2024????.nc'%(model)
    ds = xr.open_mfdataset(fp_dataset)
    stn_list = ds['station'].values.tolist()
#   Find station index
    try:
        istn = stn_list.index(station)
    except ValueError:
        print("%s not in list of: "%(station), stn_list)
        sys.exit()
    
    h    = ds['H'].values
    h    = np.squeeze(h[0,:,istn])/1000.

    Species = None
    speciestitle = ''
    if(species == 'du'):
        Species = ['DU']
        speciestitle = 'Dust '
    if(species == 'su'):
        Species = ['SU']
        speciestitle = 'Sulfate '
    if(species == 'cc'):
        Species = ['OC','BC']
        if(model == 'res'):
            Species = ['OC','BR','BC']
        speciestitle = 'Carbonaceous '

    if varn != 'EXT532nm':
        varv = ds[varn].values
        varv = np.transpose(np.squeeze(varv[:,:,istn])*scale)
        clevs = np.arange(y0,y1,(y1-y0)/50.)
        print(varv.shape,np.max(varv))
    else:
        varv = ds[varn].values
        varv = np.log10(np.transpose(varv[:,:,istn]*scale))
        clevs = np.arange(-2,0.,.02)
        title = 'Extinction'
        units = '532 nm, km-1'
        print(varv.shape)
        
    unitstr = ''
    if units != '':
        unitstr = '[%s]'%(units)

    titlestr = title
    if title == '':
        titlestr = varn

    fig, ax = plt.subplots(figsize=(20, 6))
    time  = ds.time.values
    ntime = ds.dims['time']

    ax.set_ylim(0,12)
    plt.ylabel('Altitude [km]')
    dtFmt = mdates.DateFormatter('%m/%d') # define the formatting
    plt.gca().xaxis.set_major_formatter(dtFmt) # apply the format to the desired axis
    cf  = ax.contourf(time,h, varv,clevs,cmap=cm,extend='max')
    if varn != 'EXT' and varn != 'EXT532nm':
        cbar1 = plt.colorbar(cf,ax=ax)
    else:
        cbar1 = plt.colorbar(cf,ax=ax,ticks=[-2,-1,0],
                             format=mticker.FixedFormatter(['0.01', '0.1', '1']) )
    cbar1.set_label(label='%s %s'%(titlestr,unitstr),
                    size=16,rotation=270.,labelpad=25)
    ax.set_facecolor('black')
    plt.title('%s, %s'%(station,titlestr), size=20)
    plt.savefig('stn_plots/%s.%s.%s.png'%(model,station,varn))
    plt.close(fig)

if __name__ == "__main__":
    stations = ['Resolute_Bay', 'Pituffik', 'Kaffeklubben_Island', 'Eureka',
                'Alert', 'Baffin_Bay', 'Belle', 'Fram_Strait', 'Lisee',
                'Northern_Arctic_Ocean', 'Station_Nord', 'Summit_Greenland',
                'Svalbard_Zeppelin']
    for station in stations:
        print(station)
        model   = 'MERRA2.inst3_3d_aer_Nv'
        plot(model=model,station=station,varn='T',units='K',title='Temperature',y0=200,y1=300)
        plot(model=model,station=station,varn='EXT532nm')

