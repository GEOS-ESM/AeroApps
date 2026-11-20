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
from optparse   import OptionParser   # Command-line args

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

def plotext(ictFile,model="fp",collection="inst3-3d-AER-Nv"):

#   Get the ICARTT file describing the trajectory
    if ictFile.find('P3B') > 0:
        aircraft = 'P3B'
        i0 = ictFile.find('P3B')+4
    if ictFile.find('G3')  > 0:
        aircraft = 'G3'
        i0 = ictFile.find('G3')+3
    if ictFile.find('Learjet') > 0:
        aircraft = 'Learjet'
        i0 = ictFile.find('Learjet')+8
    m = ICARTT(ictFile)
    alt, lon, lat, tyme = m.Nav['Altitude'], m.Nav['Longitude'], m.Nav['Latitude'], m.Nav['Time']
    yyyymmdd = ictFile[i0:i0+11]
    dateout  = ictFile[i0:i0+4]+"-"+ictFile[i0+4:i0+6]+"-"+ictFile[i0+6:i0+8]
    print(ictFile, aircraft, yyyymmdd)

    modeltitle='GEOS-FP'
    modname = model
    if(model == 'fp'):
        modname = "GEOS-FP"

#   Get the sampled file
    dirname = f"samples/ARCSIX/sampled/{aircraft}/{modname}/{dateout}"
    sampleFile = f"./{dirname}/ARCSIX-{modname}-{collection}-{aircraft}_Model_{yyyymmdd}.nc"
    config = '/home/pcolarco/silo/GMAOpyobs/src/config/m2_pm25.yaml'

    if(model == 'res'):
        sampleFile = '%s.inst3d_aer_v.%s.'%(model,aircraft)+yyyymmdd+'.nc'
        config = './v2xx_optics.yaml'
    print(sampleFile)
    optics = G2GAOP(sampleFile,config=config)
    
#   Get the extinction profile
    ext = optics.getAOPext(wavelength=532)

#   Get the height
    asmfile = sampleFile
    if isinstance(asmfile,xr.Dataset):
        asm = asmfile
    else:
        asm = xc.open_mfdataset(asmfile)
    z   = asm['H']

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
    ntime = ext.sizes['time']
    nlev = ext.sizes['level']
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
    
    plt.title('Aerosol Concentration, %s track: '%(aircraft)+dateout, size=20)
    ofname = f"{dirname}/ARCSIX-{modname}-{collection}-{aircraft}_Model_{yyyymmdd}.mass_along.png"
    plt.savefig(ofname)
    plt.close(fig)

if __name__ == "__main__":

    parser = OptionParser(usage="Usage: %prog [options] modelname date0",
                          version='xxx' )
    (options, args) = parser.parse_args()
 
#  GET ICT FILE FROM INPUT ARGUMENT LIST
    if len(args) == 2:
        ict        = args[0]
        model      = args[1]
    else:
        parser.error("must have 0 argument: icartt filename")

    plotext(ict,model=model)
