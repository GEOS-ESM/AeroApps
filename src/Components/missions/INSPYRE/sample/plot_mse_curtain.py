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

def plotext(yyyymmdd,species='cc',model='fp',aircraft='GV'):

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
    if(species == 'oc'):
        Species = ['OC']
        if(model == 'res'):
            Species = ['OC','BR']
        speciestitle = 'Organics '

#   Get the ICARTT file for the aircraft altitude
    ictFile = '/home/pcolarco/INSPYRE/data/INSPYRE-MetNav_%s_'%(aircraft)+yyyymmdd+'_R0.ict'
    if(aircraft == 'Learjet'):
        ictFile = '/home/pcolarco/INSPYRE/data/INSPYRE-NAVM300_Learjet_'+yyyymmdd+'_R0.ict'
    print(ictFile)
    m = ICARTT(ictFile)
    alt, lon, lat, tyme = m.Nav['Altitude'], m.Nav['Longitude'], m.Nav['Latitude'], m.Nav['Time']

#   Get the sampled file
    sampleFile = '%s.inst3_3d_aer_Nv.%s.'%(model,aircraft)+yyyymmdd+'.nc'
    config = '/home/pcolarco/silo/GMAOpyobs/src/config/m2_pm25.yaml'

    if(model == 'res'):
        sampleFile = '%s.inst3d_aer_v.%s.'%(model,aircraft)+yyyymmdd+'.nc'
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
    ext = optics.getAOPext(wavelength=550,Species=Species,fixrh=0.)

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
    if(species == 'oc'):
        aer = (optics.aer['OCPHOBIC']+optics.aer['OCPHILIC'])*optics.aer['AIRDENS']
        if(model == 'res'):
            aer = (optics.aer['OCPHILIC']+optics.aer['BRPHILIC']+
                   optics.aer['OCPHOBIC']+optics.aer['BRPHOBIC'])*optics.aer['AIRDENS']


    fig, ax = plt.subplots(figsize=(20, 6))
    time = ext.time.values
    ntime = ext.dims['time']
    nlev = ext.dims['lev']
    time = np.repeat(time.reshape(ntime,1),nlev,axis=1)
    ax.set_ylim(0,12)
    plt.ylabel('GPS Altitude [km]')
    dtFmt = mdates.DateFormatter('%H:%M') # define the formatting
    plt.gca().xaxis.set_major_formatter(dtFmt) # apply the format to the desired axis
    clevs = np.arange(0,10,.2)
    im  = ax.contourf(time,z/1000.,ext.SCA/aer.values/1e6,clevs,cmap=cm,extend='max')
    cf  = ax.plot(tyme,alt/1000.,color='magenta',linewidth=4)
    if(model != 'res'):
        im2 = ax.contour(time,z/1000.,cld,[0.49,0.5],colors='slategray')
    ax.set_facecolor('black')

    for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
             ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(20)

    cbar1 = plt.colorbar(im,ax=ax,ticks=[0,5,10],
                         format=mticker.FixedFormatter([ '0', '5', '10']) )
    cbar1.ax.tick_params(labelsize=16)
    cbar1.set_label(label='Mass Scattering Efficiency [m2 g-1]',
                    size=16,rotation=270.,labelpad=25)
    plt.title('%s %sMass Scattering Efficiency, %s track: '%(modeltitle,speciestitle,aircraft)+yyyymmdd, size=20)
    species_ = ''
    if(species != None):
        species_ = species+'_'
    plt.savefig('INSPYRE-MetNav_%s_'%(aircraft)+yyyymmdd+'_R0.%s_%smse_curtain.png'%(model,species_))
    plt.close(fig)

if __name__ == "__main__":

    plotext('20240607',aircraft='ER2')
    plotext('20240607',aircraft='ER2',species='oc')
    plotext('20240607',aircraft='ER2',species='su')
    plotext('20240607',aircraft='ER2',species='cc')
    sys.exit()
#   Learjet
    mmdd = ['0725','0729','0730','0801','0802','0807','0808']
    for date in mmdd:
        plotext('2024%s'%(date),aircraft='Learjet')
#        plotext('2024%s'%(date),model='geosit',aircraft='Learjet')
#        plotext('2024%s'%(date),model='res',aircraft='Learjet')

#   ER2
    mmdd = ['0524','0528','0530','0531','0603','0605','0606','0607','0610','0611','0613',
            '0722','0725','0729','0730','0801','0802','0807','0808','0809','0815']
    for date in mmdd:
        plotext('2024%s'%(date),aircraft='ER2')
#        plotext('2024%s'%(date),model='geosit',aircraft='ER2')
#        plotext('2024%s'%(date),model='res',aircraft='ER2')

    sys.exit()
#   GV
    mmdd = ['0530','0531','0603','0605','0606','0607','0610','0611','0613',
            '0806','0807','0808','0809','0815']
    for date in mmdd:
        plotext('2024%s'%(date))
        plotext('2024%s'%(date),model='geosit')
        plotext('2024%s'%(date),model='res')
