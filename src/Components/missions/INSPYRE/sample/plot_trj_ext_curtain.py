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
import datetime
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


# Plot the aerosol species mass concentration along a flight track
# for an already computed trajectory file
def plotext(ictFile,model="fp",collection="inst3-3d-AER-Nv",species=None):

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
    yyyymmdd = ictFile[i0:i0+11]
    dateout  = ictFile[i0:i0+4]+"-"+ictFile[i0+4:i0+6]+"-"+ictFile[i0+6:i0+8]
    print(ictFile, aircraft, yyyymmdd)

    modeltitle='GEOS-FP'
    modname = model
    if(model == 'fp'):
        modname = "GEOS-FP"

    Species = None
    speciestitle = ''
    if(species == 'du'):
        Species = ['DU']
        speciestitle = 'Dust '
    if(species == 'ss'):
        Species = ['SS']
        speciestitle = 'Sea Salt '
    if(species == 'su'):
        Species = ['SU']
        speciestitle = 'Sulfate '
    if(species == 'cc'):
        Species = ['OC','BR','BC']
        speciestitle = 'Carbonaceous '

#   Get the ICARTT file for the aircraft altitude
    m = ICARTT(ictFile)
    alt, lon, lat, tyme = m.Nav['Altitude'], m.Nav['Longitude'], m.Nav['Latitude'], m.Nav['Time']

#   Get the sampled file
    dirname = f"samples/INSPYRE/sampled/{aircraft}/{modname}/{dateout}"
    sampleFile = f"./{dirname}/INSPYRE-{modname}-{collection}-{aircraft}_Model_{yyyymmdd}.nc"
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
    
    fig, ax = plt.subplots(figsize=(20, 6))
    time = ext.time.values
    ntime = ext.sizes['time']
    nlev = ext.sizes['level']
    time = np.repeat(time.reshape(ntime,1),nlev,axis=1)
    plt.ylabel('GPS Altitude [km]')
    dtFmt = mdates.DateFormatter('%H:%M') # define the formatting
    plt.gca().xaxis.set_major_formatter(dtFmt) # apply the format to the desired axis
    clevs = np.arange(-2,0.,.02)
    im  = ax.contourf(time,z/1000.,np.log10(ext.EXT),clevs,cmap=cm,extend='max')
    cf  = ax.plot(tyme,alt/1000.,color='magenta',linewidth=4)
    y = np.squeeze(z[:,71]/1000.)
    im2 = ax.plot(tyme,y,color="grey",lw=2)
    d = np.zeros(len(tyme))
    ax.fill_between(tyme,y,where=y>=d,color="beige")
    if aircraft == "ER2":
        ax.set_ylim(0,22)
    else:
        ax.set_ylim(0,16)
#    ax.set_xlim([datetime.datetime(2024, 6, 7, 15, 45,0), datetime.datetime(2024, 6, 7, 15, 27,0)])
#    if(model != 'res'):
#        im2 = ax.contour(time,z/1000.,cld,[0.49,0.5],colors='slategray')
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
    if species == None:
        species = "Total"
    ofname = f"{dirname}/INSPYRE-{modname}-{collection}-{aircraft}_Model_{yyyymmdd}.{species}_extinction_curtain.png"
#    ofname = f"INSPYRE-{modname}-{collection}-{aircraft}_Model_{yyyymmdd}.{species}_extinction_curtain.png"
    print(ofname)
#    plt.show()
    plt.savefig(ofname)
    plt.close(fig)

if __name__ == "__main__":

    parser = OptionParser(usage="Usage: %prog [options] modelname date0",
                          version='xxx' )
    (options, args) = parser.parse_args()
 
#  GET OMI FILE FROM INPUT ARGUMENT LIST
    if len(args) == 2:
        ict        = args[0]
        model      = args[1]
    else:
        parser.error("must have 0 argument: icartt filename")
        
    plotext(ict,model=model)
    plotext(ict,model=model,species='cc')
#    plotext(ict,model=model,species='du')
#    plotext(ict,model=model,species='ss')
    plotext(ict,model=model,species='su')



    sys.exit()
