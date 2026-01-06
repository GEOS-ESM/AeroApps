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

def plotext(yyyymmdd,model='fp',aircraft='G3'):

    modeltitle='GEOS-FP'
    if(model == 'geosit'):
        modeltitle='GEOS-IT'
    if(model == 'res'):
        modeltitle='RESEARCH'

#   Get the ICARTT file for the aircraft altitude
    ictFile = '/home/pcolarco/ARCSIX/data/ARCSIX-MetNav_%s_'%(aircraft)+yyyymmdd+'_R0.ict'
    if(aircraft == 'Learjet'):
        ictFile = '/home/pcolarco/ARCSIX/data/ARCSIX-NAVM300_Learjet_'+yyyymmdd+'_R0.ict'
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
            z   = asm['H'].values
    else:
        asmfile = '%s.inst3_3d_asm_Nv.%s.'%(model,aircraft)+yyyymmdd+'.nc'
        if isinstance(asmfile,xr.Dataset):
            asm = asmfile
        else:
            asm = xc.open_mfdataset(asmfile)
            z   = asm['H'].values
    
#   Get the extinction profile
    ext = optics.getAOPext(wavelength=532)

#   Get the aerosol mass concentrations
    aer = ((optics.aer['DU001']+optics.aer['DU002']+optics.aer['DU003']+
           optics.aer['DU004']+optics.aer['DU005'])*optics.aer['AIRDENS']).values
    z   = np.flip(z,axis=1)
    aer = np.flip(aer,axis=1)
    aerdu = np.zeros(len(tyme))
    for itime in range(len(tyme)):
        aerdu[itime] = np.interp(alt[itime],z[itime,:],aer[itime,:])


    aer = ((optics.aer['SO4'])*optics.aer['AIRDENS']).values
#    z   = np.flip(z,axis=1)
    aer = np.flip(aer,axis=1)
    aersu = np.zeros(len(tyme))
    for itime in range(len(tyme)):
        aersu[itime] = np.interp(alt[itime],z[itime,:],aer[itime,:])

    aer = ((optics.aer['BCPHILIC']+optics.aer['OCPHILIC']+
           optics.aer['BCPHOBIC']+optics.aer['OCPHOBIC'])*optics.aer['AIRDENS']).values
    if(model == 'res'):
        aer = ((optics.aer['BCPHILIC']+optics.aer['OCPHILIC']+optics.aer['BRPHILIC']+
               optics.aer['BCPHOBIC']+optics.aer['OCPHOBIC']+optics.aer['BRPHOBIC'])*optics.aer['AIRDENS']).values
    aer = np.flip(aer,axis=1)
    aercc = np.zeros(len(tyme))
    for itime in range(len(tyme)):
        aercc[itime] = np.interp(alt[itime],z[itime,:],aer[itime,:])


    fig, ax = plt.subplots(figsize=(20, 6))
    time = ext.time.values
    ntime = ext.dims['time']
    nlev = ext.dims['lev']
    time = np.repeat(time.reshape(ntime,1),nlev,axis=1)
#    ax.set_ylim(0,12)
    plt.ylabel('GPS Altitude [km]')
    dtFmt = mdates.DateFormatter('%H:%M') # define the formatting
    plt.gca().xaxis.set_major_formatter(dtFmt) # apply the format to the desired axis
    cf,  = ax.plot(tyme,alt/1000.,color='magenta',linewidth=4,label='Altitude')

    for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
             ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(20)

    ax2 = ax.twinx()
    p1, = ax2.plot(tyme,aersu*1e9,color='green',label='Sulfate')
    p2, = ax2.plot(tyme,aercc*1e9,color='black',label='Carbonaceous')
    ax2.set(ylabel='Sulfate / Carbonaceous [ug m-3]')
    for item in ([ax2.title, ax2.xaxis.label, ax2.yaxis.label] +
             ax2.get_xticklabels() + ax2.get_yticklabels()):
        item.set_fontsize(16)

    ax3 = ax.twinx()
    ax3.spines.right.set_position(("axes", 1.075))
    p3, = ax3.plot(tyme,aerdu*1e9,color='orange',label='Dust')
    ax3.set(ylabel='Dust [ug m-3]')
    for item in ([ax3.title, ax3.xaxis.label, ax3.yaxis.label] +
             ax3.get_xticklabels() + ax3.get_yticklabels()):
        item.set_fontsize(16)

#    ax.yaxis.label.set_color(cf.get_color())
#    ax2.yaxis.label.set_color(p2.get_color())
#    ax3.yaxis.label.set_color(p3.get_color())
    ax.legend(handles=[cf,p1,p2,p3])
    
    plt.title('Aerosol Concentration, %s track: '%(aircraft)+yyyymmdd, size=20)
    plt.savefig('ARCSIX-MetNav_%s_'%(aircraft)+yyyymmdd+'_R0.%s_mass_along.png'%(model))
    plt.close(fig)

if __name__ == "__main__":

#    plotext('20240607',aircraft='P3B')
#    sys.exit()

#   P3B
    mmdd = ['0524','0528','0530','0531','0603','0605','0606','0607','0610','0611','0613',
            '0722','0725','0729','0730','0801','0802','0807','0808','0809','0815']
    for date in mmdd:
        plotext('2024%s'%(date),aircraft='P3B')
#        plotext('2024%s'%(date),model='geosit',aircraft='P3B')

    sys.exit()
#   G3
    mmdd = ['0530','0531','0603','0605','0606','0607','0610','0611','0613',
            '0806','0807','0808','0809','0815']
    for date in mmdd:
        plotext('2024%s'%(date))
        plotext('2024%s'%(date),model='geosit')
        plotext('2024%s'%(date),model='res')
