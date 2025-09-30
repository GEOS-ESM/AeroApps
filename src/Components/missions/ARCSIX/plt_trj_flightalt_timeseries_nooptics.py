#!/usr/bin/env python3
#import xarray as xr
from pyobs.aop import G2GAOP

import xarray as xr
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


def plot(yyyymmdd, model='c180R_arcsix.inst3d_aer_v', varn='T',aircraft='G3',
         y0=250,y1=310.,scale=1.,units='',title='',species=None):

    fp_dataset = 'trj_samples/%s.%s.%s.nc'%(model,aircraft,yyyymmdd)
    ds  = xr.open_mfdataset(fp_dataset)
    z   = ds['H'].values
    gps = ds['GPS_ALTITUDE'].values
    
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
        varv = np.squeeze(varv[:,:])*scale
        clevs = np.arange(y0,y1,(y1-y0)/50.)
    else:
        varv = ds[varn].values*scale
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

    varvi = np.zeros((ntime))
    for i in range(0,ntime):
        varvi[i] = np.interp(gps[i],z[i,:],varv[i,:])

    plt.ylabel(titlestr)
    dtFmt = mdates.DateFormatter('%H:%M') # define the formatting
    plt.gca().xaxis.set_major_formatter(dtFmt) # apply the format to the desired axis
    cf = ax.plot(time,varvi,color='magenta',linewidth=4)
    plt.title('%s (%s), %s'%(aircraft,yyyymmdd,titlestr), size=20)
    plt.savefig('trj_plots/%s.%s.%s.%s.flightalt.png'%(model,aircraft,varn,yyyymmdd))
    plt.close(fig)

if __name__ == "__main__":
    model   = 'c180R_arcsix.inst3d_aer_v'
    model   = 'fp.inst3_3d_aer_Nv'

#   Learjet
    mmdd = ['0725','0729','0730','0801','0802','0807','0808']
    for date in mmdd:
        plot("2024%s"%(date),model=model,aircraft='Learjet',varn='T',units='K',
             title='Temperature',y0=200,y1=300)
        plot("2024%s"%(date),model=model,aircraft='Learjet',varn='EXT532nm')

#   P3B
    mmdd = ['0524','0528','0530','0531','0603','0605','0606','0607','0610','0611','0613',
            '0722','0725','0729','0730','0801','0802','0807','0808','0809','0815']
    for date in mmdd:
        plot("2024%s"%(date),model=model,aircraft='P3B',varn='T',units='K',
             title='Temperature',y0=200,y1=300)
        plot("2024%s"%(date),model=model,aircraft='P3B',varn='EXT532nm')
    
#   G3
    mmdd = ['0530','0531','0603','0605','0606','0607','0610','0611','0613',
            '0806','0807','0808','0809','0815']
    for date in mmdd:
        plot("2024%s"%(date),model=model,aircraft='G3',varn='T',units='K',
             title='Temperature',y0=200,y1=300)
        plot("2024%s"%(date),model=model,aircraft='G3',varn='EXT532nm')
