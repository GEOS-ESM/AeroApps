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

def plotext(yyyymmdd,species='cc',model='fp',aircraft='G3'):

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

#   Get the ICARTT file for the aircraft altitude
    ictFile = '/home/pcolarco/ARCSIX/data/ARCSIX-MetNav_%s_'%(aircraft)+yyyymmdd+'_R0.ict'
    if(aircraft == 'Learjet'):
        ictFile = '/home/pcolarco/ARCSIX/data/ARCSIX-NAVM300_Learjet_'+yyyymmdd+'_R0.ict'
    print(ictFile)
    m = ICARTT(ictFile)
    alt, lon, lat, tyme = m.Nav['Altitude'], m.Nav['Longitude'], m.Nav['Latitude'], m.Nav['Time']

#   Get the sampled file
    sampleFile = '%s.inst3_3d_aer_Nv.%s.'%(model,aircraft)+yyyymmdd+'.nc'
    model = 'c180R_arcsix.inst3d_aer_v'
    sampleFile = 'samples/%s.%s.%s.nc'%(model,aircraft,yyyymmdd)
    config = './g2g_pm25.yaml'

    if(model == 'res'):
        sampleFile = '%s.inst3d_aer_v.%s.'%(model,aircraft)+yyyymmdd+'.nc'
        config = './v2xx_optics.yaml'
    print(sampleFile)
    optics = G2GAOP(sampleFile,config=config)

    asmfile = sampleFile
    if isinstance(asmfile,xr.Dataset):
        asm = asmfile
    else:
        asm = xc.open_mfdataset(asmfile)
    z   = asm['H']
#    cld = asm['CLOUD']
    
#   Get the extinction profile
    ext = optics.getAOPext(wavelength=532,Species=Species)

#   Get the aerosol mass concentrations
    if(species == 'du'):
        aer = (optics.aer['DU001']+optics.aer['DU002']+optics.aer['DU003']+
               optics.aer['DU004']+optics.aer['DU005'])*optics.aer['AIRDENS']
    if(species == 'su'):
        aer = (optics.aer['SO4'])*optics.aer['AIRDENS']
    if(species == 'cc'):
        aer = (optics.aer['BCPHILIC']+optics.aer['OCPHILIC']+
               optics.aer['BCPHOBIC']+optics.aer['OCPHOBIC'])*optics.aer['AIRDENS']
        if(model == 'res'):
            aer = (optics.aer['BCPHILIC']+optics.aer['OCPHILIC']+optics.aer['BRPHILIC']+
                   optics.aer['BCPHOBIC']+optics.aer['OCPHOBIC']+optics.aer['BRPHOBIC'])*optics.aer['AIRDENS']


    fig, ax = plt.subplots(figsize=(20, 6))
    time = ext.time.values
    ntime = ext.dims['time']
    nlev = ext.dims['lev']
    time = np.repeat(time.reshape(ntime,1),nlev,axis=1)
    ax.set_ylim(0,12)
    plt.ylabel('GPS Altitude [km]')
    dtFmt = mdates.DateFormatter('%H:%M') # define the formatting
    plt.gca().xaxis.set_major_formatter(dtFmt) # apply the format to the desired axis
    clevs = np.arange(-1,1,.02)
    im  = ax.contourf(time,z/1000.,np.log10(aer.values*1e9),clevs,cmap=cm,extend='max')
    print(time.shape)
    print(z.shape)
    print(aer.shape)
    cf  = ax.plot(tyme,alt/1000.,color='magenta',linewidth=4)
#    if(model != 'res'):
#        im2 = ax.contour(time,z/1000.,cld,[0.49,0.5],colors='slategray')
    ax.set_facecolor('black')

    for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
             ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(20)

    cbar1 = plt.colorbar(im,ax=ax,ticks=[-1,0,1],
                         format=mticker.FixedFormatter([ '0.1', '1', '10']) )
    cbar1.ax.tick_params(labelsize=16)
    cbar1.set_label(label='Concentration [ug m-3]',
                    size=16,rotation=270.,labelpad=25)
    plt.title('%s %sAerosol Concentration, %s track: '%(modeltitle,speciestitle,aircraft)+yyyymmdd, size=20)
    species_ = ''
    if(species != None):
        species_ = species+'_'
    plt.savefig('ARCSIX-MetNav_%s_'%(aircraft)+yyyymmdd+'_R0.%s_%smass_curtain.png'%(model,species_))
    plt.close(fig)

if __name__ == "__main__":

    plotext('20240607',aircraft='P3B')
    plotext('20240607',aircraft='P3B',species='du')
    plotext('20240607',aircraft='P3B',species='su')
    plotext('20240607',aircraft='P3B',species='cc')
    sys.exit()
#   Learjet
    mmdd = ['0725','0729','0730','0801','0802','0807','0808']
    for date in mmdd:
        plotext('2024%s'%(date),aircraft='Learjet')
#        plotext('2024%s'%(date),model='geosit',aircraft='Learjet')
#        plotext('2024%s'%(date),model='res',aircraft='Learjet')

#   P3B
    mmdd = ['0524','0528','0530','0531','0603','0605','0606','0607','0610','0611','0613',
            '0722','0725','0729','0730','0801','0802','0807','0808','0809','0815']
    for date in mmdd:
        plotext('2024%s'%(date),aircraft='P3B')
#        plotext('2024%s'%(date),model='geosit',aircraft='P3B')
#        plotext('2024%s'%(date),model='res',aircraft='P3B')

    sys.exit()
#   G3
    mmdd = ['0530','0531','0603','0605','0606','0607','0610','0611','0613',
            '0806','0807','0808','0809','0815']
    for date in mmdd:
        plotext('2024%s'%(date))
        plotext('2024%s'%(date),model='geosit')
        plotext('2024%s'%(date),model='res')
