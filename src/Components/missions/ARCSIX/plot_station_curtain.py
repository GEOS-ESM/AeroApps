#!/usr/bin/env python3
import sys
from pyobs.sampler import TRAJECTORY
from pyobs.icartt import ICARTT
from pyobs.aop import G2GAOP
from pyobs.sampler import addVertCoord

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.path as mpath
from matplotlib import cm #cm is colormap
import matplotlib.dates as mdates
import matplotlib.colors as colors
import matplotlib.ticker as mticker
import matplotlib
matplotlib.use('agg')

import xarray as xr
import pyobs.xrctl  as xc

#Get the HALO RGB colors
import csv
from matplotlib.colors import LinearSegmentedColormap
with open('/home/pcolarco/lib/halo_colorbar.csv', newline='') as csvfile:
    reader = csv.DictReader(csvfile, delimiter=',',fieldnames=['r','g','b'])
    i = 0
    rgb = []
    for row in reader:
        rgb.append([int(row['r'])/255.,int(row['g'])/255.,int(row['b'])/255.])
cm = LinearSegmentedColormap.from_list(
        'my_map', rgb)

def plotext(species=None,model='fp.inst3_3d_aer_Nv.short',station='home'):

    modeltitle='GEOS-FP'
    if(model == 'geosit'):
        modeltitle='GEOS-IT'
    if(model == 'res'):
        modeltitle='RESEARCH'

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

#   Get the sampled file
    sampleFile = './samples/%s.stations.nc'%(model)
    config = '/home/pcolarco/silo/GMAOpyobs/src/config/m2_pm25.yaml'

    if(model == 'res'):
        config = './v2xx_optics.yaml'
    print(sampleFile)
    optics = G2GAOP(sampleFile,config=config)

    if(model=='res'):
        asmfile = sampleFile
        if isinstance(asmfile,xr.Dataset):
            asm = asmfile
        else:
            asm = xc.open_mfdataset(asmfile)
            z   = asm['H']
    else:
        asmfile = '%s.inst3_3d_asm_Nv.%s.'%(model,aircraft)+yyyymmdd+'.nc'
        if isinstance(asmfile,xr.Dataset):
            asm = asmfile
        else:
            asm = xc.open_mfdataset(asmfile)
            z   = asm['H']
            cld = asm['CLOUD']
    
#   Get the extinction profile
    ext = optics.getAOPext(wavelength=532,Species=Species)
    
    fig, ax = plt.subplots(figsize=(20, 6))
    time = ext.time.values
    ntime = ext.dims['time']
    nlev = ext.dims['lev']
    time = np.repeat(time.reshape(ntime,1),nlev,axis=1)
    ax.set_ylim(0,12)
    plt.ylabel('GPS Altitude [km]')
    dtFmt = mdates.DateFormatter('%H:%M') # define the formatting
    plt.gca().xaxis.set_major_formatter(dtFmt) # apply the format to the desired axis
    clevs = np.arange(-2,0.,.02)
    im  = ax.contourf(time,z/1000.,np.log10(ext.EXT),clevs,cmap=cm,extend='max')
    cf  = ax.plot(tyme,alt/1000.,color='magenta',linewidth=4)
    if(model != 'res'):
        im2 = ax.contour(time,z/1000.,cld,[0.49,0.5],colors='slategray')
    ax.set_facecolor('black')

    for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
             ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(20)

    cbar1 = plt.colorbar(im,ax=ax,ticks=[-2,-1,0],
                         format=mticker.FixedFormatter(['0.01', '0.1', '1']) )
    cbar1.ax.tick_params(labelsize=16)
    cbar1.set_label(label='%s %sAerosol Extinction [532 nm, km-1]'%(modeltitle,speciestitle),
                    size=16,rotation=270.,labelpad=25)
    plt.title('%s track: '%(aircraft)+yyyymmdd, size=20)
    species_ = ''
    if(species != None):
        species_ = species+'_'
    plt.savefig('ARCSIX-MetNav_%s_'%(aircraft)+yyyymmdd+'_R0.%s_%sext_curtain.png'%(model,species_))
    plt.close(fig)

if __name__ == "__main__":

    plotext(station='home')
