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
import os
import xarray as xr
import pyobs.xrctl  as xc
from datetime import timedelta
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
    if ictFile.find('ER2') > 0:
        aircraft = 'ER2'
        i0 = ictFile.find('ER2')+4
    if ictFile.find('GV')  > 0:
        aircraft = 'GV'
        i0 = ictFile.find('GV')+3
    if ictFile.find('Learjet') > 0:
        aircraft = 'Learjet'
        i0 = ictFile.find('Learjet')+8
    m = ICARTT(ictFile)
    alt, lon, lat, tyme = m.Nav['Altitude'], m.Nav['Longitude'], m.Nav['Latitude'], m.Nav['Time']
    yyyymmdd = ictFile[i0:-4]
    dateout  = ictFile[i0:i0+4]+"-"+ictFile[i0+4:i0+6]+"-"+ictFile[i0+6:i0+8]
    print(ictFile, aircraft, yyyymmdd)

    modeltitle='GEOS-FP'
    modname = model
    if(model == 'fp'):
        modname = "GEOS-FP"

#   Get the sampled file
    dirname = f"samples/INSPYRE/sampled/{aircraft}/{modname}/{dateout}"
    sampleFile = f"./{dirname}/INSPYRE-{modname}-{collection}-{aircraft}_Model_{yyyymmdd}.nc"
    config = './g2g_pm25.yaml'

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
    z   = np.flip(asm['H'],axis=1)

    aer = ((optics.aer['SO4'])*optics.aer['AIRDENS']).values
    aer = np.flip(aer,axis=1)
    aersu = np.zeros(len(tyme))
    for itime in range(len(tyme)):
        aersu[itime] = np.interp(alt[itime],z[itime,:],aer[itime,:])

    aer = ((optics.aer['BCPHILIC']+optics.aer['OCPHILIC']+optics.aer['BRPHILIC']+
           optics.aer['BCPHOBIC']+optics.aer['OCPHOBIC']+optics.aer['BRPHOBIC'])*optics.aer['AIRDENS']).values
    aer = np.flip(aer,axis=1)
    aercc = np.zeros(len(tyme))
    for itime in range(len(tyme)):
        aercc[itime] = np.interp(alt[itime],z[itime,:],aer[itime,:])

    aer = ((optics.aer['BCPHILIC']+optics.aer['BCPHOBIC'])*optics.aer['AIRDENS']).values
    aer = np.flip(aer,axis=1)
    aerbc = np.zeros(len(tyme))
    for itime in range(len(tyme)):
        aerbc[itime] = np.interp(alt[itime],z[itime,:],aer[itime,:])

    fig, ax = plt.subplots(figsize=(20, 6))
    time = ext.time.values
    ntime = ext.sizes['time']
    nlev = ext.sizes['level']
    time = np.repeat(time.reshape(ntime,1),nlev,axis=1)
    if aircraft == "ER2":
        ax.set_ylim(0,22)
    else:
        ax.set_ylim(0,16)
    plt.ylabel('GPS Altitude [km]')
    dtFmt = mdates.DateFormatter('%H:%M') # define the formatting
    plt.gca().xaxis.set_major_formatter(dtFmt) # apply the format to the desired axis
    cf,  = ax.plot(tyme,alt/1000.,color='magenta',linewidth=4,label='Altitude')
    for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
             ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(20)

    ax2 = ax.twinx()
    p1, = ax2.plot(tyme,aersu*1e9,color='green',label='Sulfate',zorder=100)
    p2, = ax2.plot(tyme,aercc*1e9,color='red',label='Carbonaceous',zorder=100)
    ax2.set(ylabel='Sulfate / Carbonaceous [ug m-3]')
    for item in ([ax2.title, ax2.xaxis.label, ax2.yaxis.label] +
             ax2.get_xticklabels() + ax2.get_yticklabels()):
        item.set_fontsize(16)

    ax3 = ax.twinx()
    ax3.spines.right.set_position(("axes", 1.075))
    p3, = ax3.plot(tyme,aerbc*1e12,color='black',label='Black Carbon',zorder=100)
    ax3.set(ylabel="Black Carbon [ng m-3]")
    ax3.set_yscale("log")
    for item in ([ax3.title, ax3.xaxis.label, ax3.yaxis.label] +
             ax3.get_xticklabels() + ax3.get_yticklabels()):
        item.set_fontsize(16)



#   See if there is SP2 file to read
    if aircraft == "GV":
        fn = f"../../data/INSPYRE-SP2-BC-1Hz_GV_{yyyymmdd}.ict"
        if os.path.exists(fn):
            x = ICARTT(fn)
#           Hack for time offset in first files
            p4, = ax3.plot(x.tyme+timedelta(hours=6),x.BC_mass_90_550_nm_amb,color="dimgray",label="SP2 BC Mass",zorder=50)
#            p4, = ax3.plot(x.tyme,x.BC_mass_90_550_nm_amb,color="black",label="SP2 BC Mass")
            ax.legend(handles=[cf,p1,p2,p3,p4])
        else:
                ax.legend(handles=[cf,p1,p2,p3])
    else:
            ax.legend(handles=[cf,p1,p2,p3])

    plt.title('Aerosol Concentration, %s track: '%(aircraft)+dateout, size=20)
    ofname = f"{dirname}/INSPYRE-{modname}-{collection}-{aircraft}_Model_{yyyymmdd}.mass_along.png"
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
